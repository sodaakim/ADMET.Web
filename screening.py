from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint
from rdkit.Chem.Draw import rdMolDraw2D
import hashlib

screening_route = Blueprint('screening_route', __name__)

def custom_round(value):
    if value // 10 >= 1:
        return round(value, 2) # 정수 부분이 10 이상이면 소수점 둘째자리에서 반올림
    else:
        return round(value, 3) # 정수 부분이 10 미만이면 소수점 셋째자리에서 반올림


@screening_route.route('/screening', methods=['GET'])
def show_screening_page():
    return render_template('screening.html')


@screening_route.route('/screening', methods=['POST'])
def process_screening():
    data = request.json
    smiles_list = data.get('smiles_list', [])
    results = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            smiles_hash = hashlib.md5(smiles.encode()).hexdigest()
            filename = f"molecule_{smiles_hash}.svg"
            image_path = os.path.join('static', 'images', filename)
            drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
            # 그리기 옵션 설정
            opts = drawer.drawOptions()

            opts.atomLabelFontSize = 15
            opts.bondLineWidth = 3.0
            opts.padding = 0
            opts.multipleBondOffset = 0.15
            opts.highlightColour = (0.8, 0.8, 0.2)

            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText().replace('svg:', '')
            with open(image_path, 'w') as svg_file:
                svg_file.write(svg)
            image_url = url_for('static', filename=f'images/{filename}')

            # 분자 속성 계산
            properties = {
                'Image URL': image_url,
                'SMILES': smiles,

                # Physicochemical Property
                'Formula': rdMolDescriptors.CalcMolFormula(mol),
                'Molecular Weight': custom_round(Descriptors.MolWt(mol)),
                'Volume': custom_round(Descriptors.MolMR(mol)),
                'Density': custom_round(Descriptors.MolWt(mol) / Descriptors.MolMR(mol)) if Descriptors.MolMR(mol) else 0 ,
                'Num. heavy atoms': rdMolDescriptors.CalcNumHeavyAtoms(mol),
                'Num. arom. heavy atoms': len([atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]),
                'formal Charge': Chem.GetFormalCharge(mol),
                'TPSA': custom_round(rdMolDescriptors.CalcTPSA(mol)),

                'logP': custom_round(Descriptors.MolLogP(mol)),

                'Num. H-bond acceptors': rdMolDescriptors.CalcNumHBA(mol),
                'Num. H-bond donors': rdMolDescriptors.CalcNumHBD(mol),
                'Num. rotatable bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),

            }
            results.append(properties)
        else:
            results.append({'error': 'Invalid SMILES string', 'smiles': smiles})

    session['smiles_data'] = results
    return jsonify({'message': 'Properties calculated and saved', 'results': results})




'''
steps = [
        {'title': 'Step 1', 'description': 'Input Smiles to TextArea'},
        {'title': 'Step 2', 'description': 'Click "Send" Button before analyze'},
        {'title': 'Step 3', 'description': 'Click "Analyze Smiles" Button'}
    ]
    cardsA = [
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'}
    ]
    cardsI = [
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'}
    ]
'''