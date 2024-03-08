from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint
from rdkit.Chem.Draw import rdMolDraw2D
import hashlib

# Blueprint 객체 생성
analyze_route = Blueprint('analyze_route', __name__)

def custom_round(value):
    # 정수 부분이 10 이상인지 체크
    if value // 10 >= 1:
        # 정수 부분이 10 이상이면 소수점 둘째자리에서 반올림
        return round(value, 2)
    else:
        # 정수 부분이 10 미만이면 소수점 셋째자리에서 반올림
        return round(value, 3)

@analyze_route.route('/analyze', methods=['GET', 'POST'])
def analyze():
    if request.method == 'POST':
        smiles = request.json['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            smiles_hash = hashlib.md5(smiles.encode()).hexdigest()
            filename = f"molecule_{smiles_hash}.svg"
            image_path = os.path.join('static', 'images', filename)
            drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
            # 그리기 옵션 설정
            opts = drawer.drawOptions()

            opts.useBWAtomPalette()  # 흑백으로 원자 표시
            opts.atomLabelFontSize = 15  # 원자 라벨의 글꼴 크기
            opts.bondLineWidth = 2.0  # 본드 선의 두께
            opts.padding = 0.2  # 이미지 가장자리의 패딩
            opts.multipleBondOffset = 0.15  # 다중 결합의 오프셋
            opts.highlightColour = (0.8, 0.8, 0.2)  # 하이라이트 색상 (RGB)

            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText().replace('svg:', '')
            with open(image_path, 'w') as svg_file:
                svg_file.write(svg)
            image_url = url_for('static', filename=f'images/{filename}')


            # 분자 속성 계산 후 세션에 저장
            session['properties'] = {
                # Physicochemical Property
                'Formula': rdMolDescriptors.CalcMolFormula(mol),
                'Molecular Weight': custom_round(Descriptors.MolWt(mol)),
                'Volume': custom_round(Descriptors.MolMR(mol)),
                'Density': custom_round(Descriptors.MolWt(mol) / Descriptors.MolMR(mol)) if Descriptors.MolMR(mol) else 0,
                'Num. heavy atoms': rdMolDescriptors.CalcNumHeavyAtoms(mol),
                'Num. arom. heavy atoms': len([atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]),
                'formal Charge': Chem.GetFormalCharge(mol),
                'TPSA': custom_round(rdMolDescriptors.CalcTPSA(mol)),

                'logP': custom_round(Descriptors.MolLogP(mol)),

                'Num. H-bond acceptors': rdMolDescriptors.CalcNumHBA(mol),
                'Num. H-bond donors': rdMolDescriptors.CalcNumHBD(mol),
                'Num. rotatable bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
            }
            session['image'] = image_url
            return jsonify({'smiles': smiles, 'message': 'Properties calculated and saved'}), 200
        else:
            return jsonify({'error': 'Invalid SMILES string'}), 400
    else:
        # GET 요청 처리 로직
        steps = [
            {'title': 'Step 1', 'description': 'paste admet1'},
            {'title': 'Step 2', 'description': 'paste admet2'},
            {'title': 'Step 3', 'description': 'paste admsdsdsdsdsd'}
        ]
        return render_template('analyze.html', steps=steps)
