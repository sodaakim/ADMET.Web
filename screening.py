from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint
import hashlib

# Blueprint 객체 생성
screening_route = Blueprint('screening_route', __name__)

def custom_round(value):
    # 정수 부분이 10 이상인지 체크
    if value // 10 >= 1:
        # 정수 부분이 10 이상이면 소수점 둘째자리에서 반올림
        return round(value, 2)
    else:
        # 정수 부분이 10 미만이면 소수점 셋째자리에서 반올림
        return round(value, 3)


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
            filename = f"molecule_{smiles_hash}.png"
            image_path = os.path.join('static', 'images', filename)
            img = Draw.MolToImage(mol)
            img.save(image_path)
            image_url = url_for('static', filename=f'images/{filename}')

            properties = {
                'Formula': rdMolDescriptors.CalcMolFormula(mol),
                'Molecular Weight': custom_round(Descriptors.MolWt(mol)),
                'Image URL': image_url,
                'SMILES': smiles
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