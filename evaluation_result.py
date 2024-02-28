from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint

# Blueprint 객체 생성
evaluation_result_route = Blueprint('evaluation_result_route', __name__)

# http://127.0.0.1:8000/evaluation/%3Csmiles%3E
# smiles_input

@evaluation_result_route.route('/evaluation-result')
def evaluation_result(smiles):
    details = get_details_from_session(smiles)
    smiles = request.args.get('smiles')

    # 여러 SMILES 문자열에 대한 결과를 세션에서 가져옴
    if 'results' in session:
        results = session['results']
        # 결과를 처리하는 로직이 여기에 들어갑니다.
        # 예를 들어, 특정 SMILES 문자열에 해당하는 결과를 찾아서 사용할 수 있습니다.
        for result in results:
            if result.get('smiles') == smiles:  # 특정 SMILES 문자열에 대한 결과 찾기
                prop_values = result
                break
        else:
            prop_values = None  # 해당하는 SMILES 문자열의 결과가 없는 경우

        # 결과가 있을 경우에만 페이지에 표시할 데이터를 준비
        if prop_values:
            properties = [
                # 여기에 표시할 속성 데이터를 prop_values에서 가져와 구성합니다.
                {'title': 'Physicochemical Properties', 'contents': [
                    {'subtitle': 'Formula', 'model1': '   ', 'model2': '   ',
                     'model3': prop_values.get('Formula', 'N/A'), 'color': '#103B59',
                     'tooltipText': 'Formula: A symbolic representation...'},
                    {'subtitle': 'Molecular Weight', 'model1': '   ', 'model2': '   ',
                     'model3': prop_values.get('Molecular Weight', 'N/A'), 'color': '#103B59',
                     'tooltipText': 'Molecular Weight (MW): The mass of a molecule.'},
                    # 나머지 속성들...
                ]}
                # 추가적으로 더 많은 속성을 여기에 추가할 수 있습니다.
            ]
            return render_template('evaluation.html', properties=properties, smiles=smiles)
        else:
            # 해당하는 결과가 없는 경우, 오류 메시지를 표시할 수 있습니다.
            return "No results found for the given SMILES string."
    else:
        # 세션에 결과가 전혀 없는 경우
        return "No analysis results available."

def get_details_from_session(smiles):
    # 세션에서 특정 SMILES 문자열에 대한 세부 정보를 추출
    pass