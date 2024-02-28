from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint

# Blueprint 객체 생성
result_route = Blueprint('result_route', __name__)


def create_cards(analysis_results):
    # 세션에서 분석 결과를 가져옵니다.
    if 'results' in session:
        results = session['results']
        cards = []
        for result in results:
            # 각 결과에 대한 카드를 생성합니다.
            # 예시에서는 'image_url'과 'smiles' 키가 각 결과에 포함되어 있다고 가정합니다.
            card = {
                'img': result.get('image_url', 'default_image_url'),  # 결과에 이미지 URL이 없는 경우 대체 이미지 URL 사용
                'smiles': result.get('smiles', 'N/A')  # 결과에 SMILES 문자열이 없는 경우 'N/A' 사용
            }
            cards.append(card)
        return cards
    else:
        # 세션에 결과가 없는 경우 빈 리스트 반환
        return []

@result_route.route('/result')
def result():
    analysis_results = session.get('analysis_results', [])

    cards = create_cards(analysis_results)
    per_page = 18
    page = request.args.get('page', 1, type=int)
    total_pages = (len(cards) + per_page - 1) // per_page

    # 현재 페이지에 표시할 카드의 서브셋을 계산합니다.
    start_index = (page - 1) * per_page
    end_index = start_index + per_page
    page_cards = cards[start_index:end_index]

    cardColors = [
        '#778899', '#6892b2', '#8f8ec5', '#366d98', '#a38893', '#88a388',
        '#800020', '#4a6862', '#778899', '#88a398', '#808000', '#c0c080'
    ]

    properties = [  # property
        # content : Physicochemical Property
        {'title': 'Physicochemical Properties', 'color': cardColors[0],
         'contents': [{'subtitle': 'Formula'},          {'subtitle': 'Molecular Weight'},   {'subtitle': 'Volume'},
                      {'subtitle': 'Density'},          {'subtitle': 'Num. heavy atoms'},   {'subtitle': 'Num. arom. heavy atoms'},
                      {'subtitle': 'formal Charge'},    {'subtitle': 'TPSA'},               {'subtitle': 'logS'},
                      {'subtitle': 'logP'},             {'subtitle': 'logD'},               {'subtitle': 'DMSO'},
                      {'subtitle': 'nHA'},              {'subtitle': 'nHD'},                {'subtitle': 'nRot'},
                      {'subtitle': 'nRing'},            {'subtitle': 'nHet'},               {'subtitle': 'nRig'},
                      {'subtitle': 'pKa'} ]},

        # content : Medicinal Chemistry
        {'title': 'Medicinal Chemistry', 'color': cardColors[1],
         'contents': [{'subtitle': 'SAscore'},                  {'subtitle': 'NPscore'},    {'subtitle': 'Lipinski Rule'},
                      {'subtitle': 'Pfizer Rule'},              {'subtitle': 'Ghose'},      {'subtitle': 'GSK Rule'},
                      {'subtitle': 'Veber (GSK) filter'},    {'subtitle': 'Golden Triangle'}, {'subtitle': 'PAINS'},
                      {'subtitle': 'Leadlikeness'},             {'subtitle': 'Brenk'} ]},

        # content : Distribution
        {'title': 'Distribution', 'color': cardColors[2],
         'contents': [{'subtitle': 'PPB'},              {'subtitle': 'Fu'},     {'subtitle': 'Vd'},
                      {'subtitle': 'BBB Penetration'},  {'subtitle': 'Log Kp'}, {'subtitle': 'Organic anion'},
                      {'subtitle': 'HPP Binding'},      {'subtitle': 'GI absorption '}, {'subtitle': 'P-gp substrate '} ]},

        # content : Druglikeness
        {'title': 'Druglikeness', 'color': cardColors[3],
         'contents': [{'subtitle': 'Egan'},   {'subtitle': 'Bioavailability Score'} ]},

        # content : Absorption
        {'title': 'Absorption', 'color': cardColors[4],
         'contents': [{'subtitle': 'Caco-2 Permeability'},  {'subtitle': 'MDCK Permeability'},  {'subtitle': 'Pgp-inhibitor'},
                      {'subtitle': 'Pgp-substrate', },      {'subtitle': 'HIA'} ]},

        # content : Excretion
        {'title': 'Excretion', 'color': cardColors[5],
         'contents': [{'subtitle': 'Clearance'},    {'subtitle': 'Half-life'} ]},

        # content : Metabolism
        {'title': 'Metabolism', 'color': cardColors[6],
         'contents': [{'subtitle': 'CYP1A2 inhibitor'},     {'subtitle': 'CYP1A2 substrate'},   {'subtitle': 'CYP2C19 inhibitor'},
                      {'subtitle': 'CYP2C19 substrate'},    {'subtitle': 'CYP2C9 inhibitor'},   {'subtitle': 'CYP2C9 substrate'},
                      {'subtitle': 'CYP2D6 inhibitor'},     {'subtitle': 'CYP2D6 substrate'},   {'subtitle': 'CYP3A4 inhibitor'},
                      {'subtitle': 'CYP3A4 substrate'} ]},

        # content : Toxicity
        {'title': 'Toxicity', 'color': cardColors[7],
         'contents': [{'subtitle': 'hERG Blockers'},        {'subtitle': 'DILI'},           {'subtitle': 'AMES Toxicity'},
                      {'subtitle': 'Carcinogencity'},       {'subtitle': 'Eye Irritation'}, {'subtitle': 'Eye Corrosion'},
                      {'subtitle': 'Hepatotoxicity'},       {'subtitle': 'Rat Oral Acute Toxicity'},    {'subtitle': 'Skin Sensitization'} ]},

        # content : Cytotoxicity
        {'title': 'Cytotoxicity', 'color': cardColors[8],
         'contents': [{'subtitle': 'HepG2'},    {'subtitle': 'NIH 3T3'} ]},

        # content : Environmental Toxicity
        {'title': 'Environmental Toxicity', 'color': cardColors[9],
         'contents': [{'subtitle': 'Bioconcentration Factors'}, {'subtitle': 'Biodegration factors'},   {'subtitle': 'LD50'} ]},

        # content : Biological Pathway Toxicity
        {'title': 'Biological Pathway Toxicity', 'color': cardColors[10],
         'contents': [
             {'subtitle': 'NR-AR (Androgen receptor)'},
             {'subtitle': 'NR-AR-LBD (ligand -binding)'},
             {'subtitle': 'NR-AhR (Aryl hydrocarbon)'},
             {'subtitle': 'NR-Aromatase'},
             {'subtitle': 'NR-ER (Estrogen receptor)'},
             {'subtitle': 'NR-ER-LBD'},
             {'subtitle': 'NR-PPAR-gamma (Peroxisome)'},
             {'subtitle': 'SR-ARE (Antioxidant response)'},
             {'subtitle': 'SR-ATAD5 (ATPase family domain)'},
             {'subtitle': 'SR-HSE (Heat shock factor)'},
             {'subtitle': 'SR-MMP (Mitochondrial membrane)'},
             {'subtitle': 'SR-p53'}
         ]},
        # content : Toxicophore Rules
        {'title': 'Toxicophore Rules', 'color': cardColors[11],
         'contents': [
             {'subtitle': 'Acute Toxicity Rule'},
             {'subtitle': 'Genotoxic Carcinogenicity Rule'},
             {'subtitle': 'NonGenotoxic Carcinogenicity Rule'},
             {'subtitle': 'Skin Sensitization Rule'},
             {'subtitle': 'Aquatic Toxicity Rule'},
             {'subtitle': 'NonBiodegradable Rule'},
             {'subtitle': 'SureChEMBL Rule'},
             {'subtitle': 'FAF-Drugs4 Rule'}
         ]},
    ]
    return render_template('result.html', origin_cards=cards, cards=page_cards, page=page, total_pages=total_pages, properties=properties)