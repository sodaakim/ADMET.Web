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
                {'subtitle': 'Formula',                 'num': prop_values.get('Formula', 'N/A'),                   'color': '#103B59'},
                {'subtitle': 'Molecular Weight',        'num': prop_values.get('Molecular Weight', 'N/A'),          'color': '#103B59'},
                {'subtitle': 'Volume',                  'num': prop_values.get('Volume', 'N/A'),                    'color': '#103B59'},
                {'subtitle': 'Density',                 'num': prop_values.get('Density', 'N/A'),                   'color': '#103B59'},
                {'subtitle': 'Num. heavy atoms',        'num': prop_values.get('Num. heavy atoms', 'N/A'),          'color': '#103B59'},
                {'subtitle': 'Num. arom. heavy atoms',  'num': prop_values.get('Num. arom. heavy atoms', 'N/A'),    'color': '#103B59'},
                {'subtitle': 'formal Charge',           'num': prop_values.get('formal Charge', 'N/A'),             'color': '#103B59'},
                {'subtitle': 'TPSA',                    'num': 'num', 'color': '#103B59'},
                {'subtitle': 'logS',                    'num': 'num', 'color': '#103B59'},
                {'subtitle': 'logP',                    'num': 'num', 'color': '#103B59'},
                {'subtitle': 'logD',                    'num': 'num', 'color': '#103B59'},
                {'subtitle': 'DMSO',                    'num': 'num', 'color': '#103B59'},
                {'subtitle': 'nHA',                     'num': 'num', 'color': '#103B59'},
                {'subtitle': 'nHD',                     'num': 'num', 'color': '#103B59'},
                {'subtitle': 'nRot',                    'num': 'num', 'color': '#103B59'},
                {'subtitle': 'nRing',                   'num': 'num', 'color': '#103B59'},
                {'subtitle': 'nHet',                    'num': 'num', 'color': '#103B59'},
                {'subtitle': 'nRig',                    'num': 'num', 'color': '#103B59'},
                {'subtitle': 'pKa',                     'num': 'num', 'color': '#103B59'}
            ]},
            # content : Medicinal Chemistry
            {'title': 'Medicinal Chemistry',
             'contents': [
                {'subtitle': 'SAscore',                 'num': 'num', 'color': '#103B59'},
                {'subtitle': 'NPscore',                 'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Lipinski Rule',           'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Pfizer Rule',             'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Ghose',                   'num': 'num', 'color': '#103B59'},
                {'subtitle': 'GSK Rule',                'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Veber (GSK) filter',   'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Golden Triangle',         'num': 'num', 'color': '#103B59'},
                {'subtitle': 'PAINS',                   'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Leadlikeness',            'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Brenk',                   'num': 'num', 'color': '#103B59'}
            ]},
            # content : Distribution
            {'title': 'Distribution',
             'contents': [
                {'subtitle': 'PPB',                     'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Fu',                      'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Vd',                      'num': 'num', 'color': '#103B59'},
                {'subtitle': 'BBB Penetration',         'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Log Kp',                  'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Organic anion',           'num': 'num', 'color': '#103B59'},
                {'subtitle': 'HPP Binding',             'num': 'num', 'color': '#103B59'},
                {'subtitle': 'GI absorption ',       'num': 'num', 'color': '#103B59'},
                {'subtitle': 'P-gp substrate ',      'num': 'num', 'color': '#103B59'}
            ]},
            # content : Druglikeness
            {'title': 'Druglikeness',
             'contents': [
                {'subtitle': 'Egan',                    'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Bioavailability Score',   'num': 'num', 'color': '#103B59'}
            ]},
            # content : Absorption
            {'title': 'Absorption',
             'contents': [
                {'subtitle': 'Caco-2 Permeability',     'num': 'num', 'color': '#103B59'},
                {'subtitle': 'MDCK Permeability',       'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Pgp-inhibitor',           'num': 'num', 'color': '#103B59'},
                {'subtitle': 'Pgp-substrate',           'num': 'num', 'color': '#103B59'},
                {'subtitle': 'HIA',                     'num': 'num', 'color': '#103B59'}
            ]},
            # content : Excretion
            {'title': 'Excretion',
             'contents': [
                 {'subtitle': 'Clearance',              'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Half-life',              'num': 'num', 'color': '#103B59'}
             ]},
            # content : Metabolism
            {'title': 'Metabolism',
             'contents': [
                 {'subtitle': 'CYP1A2 inhibitor',       'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'CYP1A2 substrate',       'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'CYP2C19 inhibitor',      'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'CYP2C19 substrate',      'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'CYP2C9 inhibitor',       'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'CYP2C9 substrate',       'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'CYP2D6 inhibitor',       'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'CYP2D6 substrate',       'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'CYP3A4 inhibitor',       'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'CYP3A4 substrate',       'num': 'num', 'color': '#103B59'}
             ]},
            # content : Toxicity
            {'title': 'Biological Target-Based Toxicity',
             'contents': [
                 {'subtitle': 'hERG Blockers',          'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'DILI',                   'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'AMES Toxicity',          'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Carcinogencity',         'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Eye Irritation',         'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Eye Corrosion',          'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Hepatotoxicity',         'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Rat Oral Acute Toxicity', 'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Skin Sensitization',     'num': 'num', 'color': '#103B59'}
             ]},
            # content : Cytotoxicity
            {'title': 'Cytotoxicity',
             'contents': [
                 {'subtitle': 'HepG2',                  'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'NIH 3T3',                'num': 'num', 'color': '#103B59'}
             ]},
            # content : Environmental Toxicity
            {'title': 'Environmental Toxicity',
             'contents': [
                 {'subtitle': 'Bioconcentration Factors',   'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Biodegration factors',       'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'LD50',                       'num': 'num', 'color': '#103B59'}
             ]},
            # content : Biological Pathway Toxicity
            {'title': 'Biological Pathway Toxicity',
             'contents': [
                 {'subtitle': 'NR-AR (Androgen receptor)',          'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'NR-AR-LBD (ligand -binding)',        'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'NR-AhR (Aryl hydrocarbon)',          'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'NR-Aromatase',                       'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'NR-ER (Estrogen receptor)',          'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'NR-ER-LBD',                          'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'NR-PPAR-gamma (Peroxisome)',         'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'SR-ARE (Antioxidant response)',      'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'SR-ATAD5 (ATPase family domain)',    'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'SR-HSE (Heat shock factor)',         'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'SR-MMP (Mitochondrial membrane)',    'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'SR-p53',                             'num': 'num', 'color': '#103B59'}
             ]},
            # content : Toxicophore Rules
            {'title': 'Toxicophore Rules',
             'contents': [
                 {'subtitle': 'Acute Toxicity Rule',                'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Genotoxic Carcinogenicity Rule',     'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'NonGenotoxic Carcinogenicity Rule',  'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Skin Sensitization Rule',            'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Aquatic Toxicity Rule',              'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'NonBiodegradable Rule',              'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'SureChEMBL Rule',                    'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'FAF-Drugs4 Rule',                    'num': 'num', 'color': '#103B59'}
             ]},
        ]

        # 쿼리 파라미터에서 이미지 URL과 SMILES 문자열을 추출
        image_url = request.args.get('image_url')
        smiles = request.args.get('smiles')

        return render_template('evaluation.html', properties=properties, image_url=image_url, smiles=smiles)
    else:
        return 'Molecule not found', 404
