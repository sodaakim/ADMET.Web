from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint
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
            filename = f"molecule_{smiles_hash}.png"
            image_path = os.path.join('static', 'images', filename)
            img = Draw.MolToImage(mol)
            img.save(image_path)
            image_url = url_for('static', filename=f'images/{filename}')

            # 분자 속성 계산
            # Physicochemical Property
            mol_formula = rdMolDescriptors.CalcMolFormula(mol)
            mol_weight = custom_round(Descriptors.MolWt(mol))
            mol_volume = custom_round(Descriptors.MolMR(mol))  # Volume의 대략적 추정
            density = custom_round(mol_weight / mol_volume) if mol_volume else 0  # Density 계산
            num_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)
            num_aromatic_heavy_atoms = len([atom for atom in mol.GetAtoms() if atom.GetIsAromatic()])
            formal_charge = Chem.GetFormalCharge(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            logp = custom_round(Descriptors.MolLogP(mol))
            num_h_acceptors = rdMolDescriptors.CalcNumHBA(mol)
            num_h_donors = rdMolDescriptors.CalcNumHBD(mol)
            num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

            # 세션에 속성 저장
            session['properties'] = {
                'Formula': mol_formula,
                'Molecular Weight': mol_weight,
                'Volume': mol_volume,
                'Density': density,
                'Num. heavy atoms': num_heavy_atoms,
                'Num. arom. heavy atoms': num_aromatic_heavy_atoms,
                'formal Charge': formal_charge,
                'TPSA': custom_round(tpsa),

                'logP': custom_round(logp),

                'Num. H-bond acceptors': num_h_acceptors,
                'Num. H-bond donors': num_h_donors,
                'Num. rotatable bonds': num_rotatable_bonds
            }

            return jsonify({'image_url': image_url, 'message': 'Properties calculated and saved'}), 200
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
