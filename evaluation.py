from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint

# Blueprint 객체 생성
evaluation_route = Blueprint('evaluation_route', __name__)

# http://127.0.0.1:8000/evaluation/%3Csmiles%3E
# smiles_input

@evaluation_route.route('/evaluation')
def evaluation():
    smiles = request.args.get('smiles')
    if 'properties' in session:
        prop_values = session['properties']
        properties = [ # property
            # content : Physicochemical Property
            {'title': 'Physicochemical Properties',
             'contents': [
                {'subtitle': 'Formula',                 'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Formula', 'N/A'),                   'color': '#103B59'},
                {'subtitle': 'Molecular Weight',        'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Molecular Weight', 'N/A'),          'color': '#103B59'},
                {'subtitle': 'Volume',                  'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Volume', 'N/A'),                    'color': '#103B59'},
                {'subtitle': 'Density',                 'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Density', 'N/A'),                   'color': '#103B59'},
                {'subtitle': 'Num. heavy atoms',        'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Num. heavy atoms', 'N/A'),          'color': '#103B59'},
                {'subtitle': 'Num. arom. heavy atoms',  'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Num. arom. heavy atoms', 'N/A'),    'color': '#103B59'},
                {'subtitle': 'formal Charge',           'model1': '   ',    'model2': '   ',    'model3': prop_values.get('formal Charge', 'N/A'),             'color': '#103B59'},
                {'subtitle': 'TPSA',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'logS',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'logP',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'logD',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'DMSO',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'nHA',                     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'nHD',                     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'nRot',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'nRing',                   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'nHet',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'nRig',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'pKa',                     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
            ]},
            # content : Medicinal Chemistry
            {'title': 'Medicinal Chemistry',
             'contents': [
                {'subtitle': 'SAscore',                 'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'NPscore',                 'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Lipinski Rule',           'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Pfizer Rule',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Ghose',                   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'GSK Rule',                'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Veber (GSK) filter',   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Golden Triangle',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'PAINS',                   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Leadlikeness',            'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Brenk',                   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
            ]},
            # content : Distribution
            {'title': 'Distribution',
             'contents': [
                {'subtitle': 'PPB',                     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Fu',                      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Vd',                      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'BBB Penetration',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Log Kp',                  'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Organic anion',           'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'HPP Binding',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'GI absorption ',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'P-gp substrate ',      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
            ]},
            # content : Druglikeness
            {'title': 'Druglikeness',
             'contents': [
                {'subtitle': 'Egan',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Bioavailability Score',   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
            ]},
            # content : Absorption
            {'title': 'Absorption',
             'contents': [
                {'subtitle': 'Caco-2 Permeability',     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'MDCK Permeability',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Pgp-inhibitor',           'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'Pgp-substrate',           'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                {'subtitle': 'HIA',                     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
            ]},
            # content : Excretion
            {'title': 'Excretion',
             'contents': [
                 {'subtitle': 'Clearance',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Half-life',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
             ]},
            # content : Metabolism
            {'title': 'Metabolism',
             'contents': [
                 {'subtitle': 'CYP1A2 inhibitor',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'CYP1A2 substrate',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'CYP2C19 inhibitor',      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'CYP2C19 substrate',      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'CYP2C9 inhibitor',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'CYP2C9 substrate',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'CYP2D6 inhibitor',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'CYP2D6 substrate',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'CYP3A4 inhibitor',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'CYP3A4 substrate',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
             ]},
            # content : Toxicity
            {'title': 'Biological Target-Based Toxicity',
             'contents': [
                 {'subtitle': 'hERG Blockers',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'DILI',                       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'AMES Toxicity',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Carcinogencity',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Eye Irritation',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Eye Corrosion',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Hepatotoxicity',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Rat Oral Acute Toxicity',    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Skin Sensitization',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
             ]},
            # content : Cytotoxicity
            {'title': 'Cytotoxicity',
             'contents': [
                 {'subtitle': 'HepG2',                      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'NIH 3T3',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
             ]},
            # content : Environmental Toxicity
            {'title': 'Environmental Toxicity',
             'contents': [
                 {'subtitle': 'Bioconcentration Factors',   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Biodegration factors',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'LD50',                       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
             ]},
            # content : Biological Pathway Toxicity
            {'title': 'Biological Pathway Toxicity',
             'contents': [
                 {'subtitle': 'NR-AR (Androgen receptor)',          'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'NR-AR-LBD (ligand -binding)',        'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'NR-AhR (Aryl hydrocarbon)',          'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'NR-Aromatase',                       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'NR-ER (Estrogen receptor)',          'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'NR-ER-LBD',                          'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'NR-PPAR-gamma (Peroxisome)',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'SR-ARE (Antioxidant response)',      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'SR-ATAD5 (ATPase family domain)',    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'SR-HSE (Heat shock factor)',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'SR-MMP (Mitochondrial membrane)',    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'SR-p53',                             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
             ]},
            # content : Toxicophore Rules
            {'title': 'Toxicophore Rules',
             'contents': [
                 {'subtitle': 'Acute Toxicity Rule',                'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Genotoxic Carcinogenicity Rule',     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'NonGenotoxic Carcinogenicity Rule',  'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Skin Sensitization Rule',            'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'Aquatic Toxicity Rule',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'NonBiodegradable Rule',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'SureChEMBL Rule',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'},
                 {'subtitle': 'FAF-Drugs4 Rule',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59'}
             ]},
        ]

        # 쿼리 파라미터에서 이미지 URL과 SMILES 문자열을 추출
        image_url = request.args.get('image_url')
        smiles = request.args.get('smiles')

        return render_template('evaluation.html', properties=properties, image_url=image_url, smiles=smiles)
    else:
        return 'Molecule not found', 404
