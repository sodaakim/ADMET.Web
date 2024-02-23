from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint

# Blueprint 객체 생성
result_route = Blueprint('result_route', __name__)

@result_route.route('/result')
def result():
    cards = [
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'},
        {'img': 'url', 'smiles': 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'}
    ]

    cardColors = [
        '#6B8E23',
        '#E6E6FA',
        '#808080',
        '#C2B280',
        '#228B22',
        '#FF6347',
        '#87CEEB',
        '#D2691E',
        '#FF8C00',
        '#778899',
        '#808000',
        '#800020',
        '#000080',
        '#FF6F61'
    ]

    properties = [ # property
        # content : Physicochemical Property
        {'title': 'Physicochemical Properties', 'color': cardColors[0],
         'contents': [
            {'subtitle': 'Formula'},
            {'subtitle': 'Molecular Weight'},
            {'subtitle': 'Volume'},
            {'subtitle': 'Density'},
            {'subtitle': 'Num. heavy atoms'},
            {'subtitle': 'Num. arom. heavy atoms'},
            {'subtitle': 'formal Charge'},
            {'subtitle': 'TPSA'},
            {'subtitle': 'logS'},
            {'subtitle': 'logP'},
            {'subtitle': 'logD'},
            {'subtitle': 'DMSO'},
            {'subtitle': 'nHA'},
            {'subtitle': 'nHD'},
            {'subtitle': 'nRot'},
            {'subtitle': 'nRing'},
            {'subtitle': 'nHet'},
            {'subtitle': 'nRig'},
            {'subtitle': 'pKa'}
        ]},
        # content : Medicinal Chemistry
        {'title': 'Medicinal Chemistry', 'color': cardColors[1],
         'contents': [
            {'subtitle': 'SAscore'},
            {'subtitle': 'NPscore'},
            {'subtitle': 'Lipinski Rule'},
            {'subtitle': 'Pfizer Rule'},
            {'subtitle': 'Ghose'},
            {'subtitle': 'GSK Rule'},
            {'subtitle': 'Veber (GSK) filter'},
            {'subtitle': 'Golden Triangle'},
            {'subtitle': 'PAINS'},
            {'subtitle': 'Leadlikeness'},
            {'subtitle': 'Brenk'}
        ]},
        # content : Distribution
        {'title': 'Distribution', 'color': cardColors[2],
         'contents': [
            {'subtitle': 'PPB'},
            {'subtitle': 'Fu'},
            {'subtitle': 'Vd'},
            {'subtitle': 'BBB Penetration'},
            {'subtitle': 'Log Kp'},
            {'subtitle': 'Organic anion'},
            {'subtitle': 'HPP Binding'},
            {'subtitle': 'GI absorption '},
            {'subtitle': 'P-gp substrate '}
        ]},
        # content : Druglikeness
        {'title': 'Druglikeness', 'color': cardColors[3],
         'contents': [
            {'subtitle': 'Egan'},
            {'subtitle': 'Bioavailability Score'}
        ]},
        # content : Absorption
        {'title': 'Absorption', 'color': cardColors[4],
         'contents': [
            {'subtitle': 'Caco-2 Permeability'},
            {'subtitle': 'MDCK Permeability'},
            {'subtitle': 'Pgp-inhibitor'},
            {'subtitle': 'Pgp-substrate',},
            {'subtitle': 'HIA'}
        ]},
        # content : Excretion
        {'title': 'Excretion', 'color': cardColors[5],
         'contents': [
             {'subtitle': 'Clearance'},
             {'subtitle': 'Half-life'}
         ]},
        # content : Metabolism
        {'title': 'Metabolism', 'color': cardColors[6],
         'contents': [
             {'subtitle': 'CYP1A2 inhibitor'},
             {'subtitle': 'CYP1A2 substrate'},
             {'subtitle': 'CYP2C19 inhibitor'},
             {'subtitle': 'CYP2C19 substrate'},
             {'subtitle': 'CYP2C9 inhibitor'},
             {'subtitle': 'CYP2C9 substrate'},
             {'subtitle': 'CYP2D6 inhibitor'},
             {'subtitle': 'CYP2D6 substrate'},
             {'subtitle': 'CYP3A4 inhibitor'},
             {'subtitle': 'CYP3A4 substrate'}
         ]},
        # content : Toxicity
        {'title': 'Toxicity', 'color': cardColors[7],
         'contents': [
             {'subtitle': 'hERG Blockers'},
             {'subtitle': 'DILI'},
             {'subtitle': 'AMES Toxicity'},
             {'subtitle': 'Carcinogencity'},
             {'subtitle': 'Eye Irritation'},
             {'subtitle': 'Eye Corrosion'},
             {'subtitle': 'Hepatotoxicity'},
             {'subtitle': 'Rat Oral Acute Toxicity'},
             {'subtitle': 'Skin Sensitization'}
         ]},
        # content : Cytotoxicity
        {'title': 'Cytotoxicity', 'color': cardColors[9],
         'contents': [
             {'subtitle': 'HepG2'},
             {'subtitle': 'NIH 3T3'}
         ]},
        # content : Environmental Toxicity
        {'title': 'Environmental Toxicity', 'color': cardColors[10],
         'contents': [
             {'subtitle': 'Bioconcentration Factors'},
             {'subtitle': 'Biodegration factors'},
             {'subtitle': 'LD50'}
         ]},
        # content : Biological Pathway Toxicity
        {'title': 'Biological Pathway Toxicity', 'color': cardColors[11],
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
        {'title': 'Toxicophore Rules', 'color': cardColors[12],
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
    return render_template('result.html', properties=properties, cards=cards)
