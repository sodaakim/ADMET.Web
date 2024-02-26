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

@evaluation_result_route.route('/evaluation_result')
def evaluation():
    smiles = request.args.get('smiles')
    if 'properties' in session:
        prop_values = session['properties']
        properties = [ # property
            # content : Physicochemical Property
            {'title': 'Physicochemical Properties',
             'contents': [
                {'subtitle': 'Formula',                 'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Formula', 'N/A'),                   'color': '#103B59', 'tooltipText': 'Formula: A symbolic representation of the chemical composition of a substance using element symbols and subscripts.'},
                {'subtitle': 'Molecular Weight',        'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Molecular Weight', 'N/A'),          'color': '#103B59', 'tooltipText': 'Molecular Weight (MW): The mass of a molecule.'},
                {'subtitle': 'Volume',                  'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Volume', 'N/A'),                    'color': '#103B59', 'tooltipText': 'Volume (Van der Waals volume): The volume occupied by a molecule based on Van der Waals interactions.'},
                {'subtitle': 'Density',                 'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Density', 'N/A'),                   'color': '#103B59', 'tooltipText': 'Density = MW / Volume'},
                {'subtitle': 'Num. heavy atoms',        'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Num. heavy atoms', 'N/A'),          'color': '#103B59', 'tooltipText': 'Num. heavy atoms: The number of non-hydrogen atoms in a molecule.'},
                {'subtitle': 'Num. arom. heavy atoms',  'model1': '   ',    'model2': '   ',    'model3': prop_values.get('Num. arom. heavy atoms', 'N/A'),    'color': '#103B59', 'tooltipText': 'Num. arom. heavy atoms: The number of aromatic non-hydrogen atoms in a molecule.'},
                {'subtitle': 'formal Charge',           'model1': '   ',    'model2': '   ',    'model3': prop_values.get('formal Charge', 'N/A'),             'color': '#103B59', 'tooltipText': 'Formal Charge: The charge assigned to an atom in a molecule according to a set of rules.'},
                {'subtitle': 'TPSA',                    'model1': '   ',    'model2': '   ',    'model3': prop_values.get('TPSA', 'N/A'),                      'color': '#103B59', 'tooltipText': 'TPSA (Topological Polar Surface Area): The surface area of a molecule that is polar.'},
                {'subtitle': 'logS',                    'model1': 'num',    'model2': 'num',    'model3': 'num',                                               'color': '#103B59', 'tooltipText': 'logS (Water Solubility): The logarithm of the water solubility of a substance.'},
                {'subtitle': 'logP',                    'model1': 'num',    'model2': 'num',    'model3': prop_values.get('logP', 'N/A'),                      'color': '#103B59', 'tooltipText': 'logP (Lipophilicity): The logarithm of the partition coefficient between n-octanol and water.'},
                {'subtitle': 'logD',                    'model1': 'num',    'model2': 'num',    'model3': 'num',                                               'color': '#103B59', 'tooltipText': 'logD: The logarithm of the distribution coefficient.'},
                {'subtitle': 'DMSO',                    'model1': 'num',    'model2': 'num',    'model3': 'num',                                               'color': '#103B59', 'tooltipText': 'DMSO (Dimethyl sulfoxide): A polar aprotic solvent commonly used in organic chemistry.'},
                {'subtitle': 'nHA',                     'model1': 'num',    'model2': 'num',    'model3': prop_values.get('Num. H-bond acceptors', 'N/A'),     'color': '#103B59', 'tooltipText': 'nHA (Num. H-acceptors): The number of hydrogen bond acceptor atoms in a molecule.'},
                {'subtitle': 'nHD',                     'model1': 'num',    'model2': 'num',    'model3': prop_values.get('Num. H-bond donors', 'N/A'),        'color': '#103B59', 'tooltipText': 'nHD (Num. H-donors): The number of hydrogen bond donor atoms in a molecule.'},
                {'subtitle': 'nRot',                    'model1': 'num',    'model2': 'num',    'model3': prop_values.get('Num. rotatable bonds', 'N/A'),      'color': '#103B59', 'tooltipText': 'nRot (Num. rotatable bonds): The number of bonds in a molecule that can rotate freely.'},
                {'subtitle': 'nRing',                   'model1': 'num',    'model2': 'num',    'model3': 'num',                                               'color': '#103B59', 'tooltipText': 'nRing (Num. rings): The number of rings in a molecule.'},
                {'subtitle': 'nHet',                    'model1': 'num',    'model2': 'num',    'model3': 'num',                                               'color': '#103B59', 'tooltipText': 'nHet (Num. heteroatoms): The number of heteroatoms (non-carbon atoms) in a molecule.'},
                {'subtitle': 'nRig',                    'model1': 'num',    'model2': 'num',    'model3': 'num',                                               'color': '#103B59', 'tooltipText': 'nRig (Num. rigid bonds): The number of bonds in a molecule that are rigid and cannot rotate freely.'},
                {'subtitle': 'pKa',                     'model1': 'num',    'model2': 'num',    'model3': 'num',                                               'color': '#103B59', 'tooltipText': 'pKa: The negative logarithm of the acid dissociation constant (Ka) of a substance, representing the strength of an acid in solution.'}
            ]},
            # content : Medicinal Chemistry
            {'title': 'Medicinal Chemistry',
             'contents': [
                {'subtitle': 'SAscore',                 'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'SAscore: A scoring system used in drug discovery to assess the synthetic accessibility of a molecule.'},
                {'subtitle': 'NPscore',                 'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'NPscore: A scoring system used to evaluate the natural product-likeness of a molecule.'},
                {'subtitle': 'Lipinski Rule',           'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'Lipinski Rule\nMax of 5 H-bond donors\nMax of 5 H-bond donors\nMax of 10 H-bond acceptors\nMax of 10 H-bond acceptors\nMolecular weight<500 daltons\nMolecular weight<500 daltons\nLog P (partition coefficient)≤5\nLog P (partition coefficient)≤5'},
                {'subtitle': 'Pfizer Rule',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'Pfizer Rule:\n300≤Molecular weight≤500 daltons\n−0.4≤Log P≤5.6\n2≤H-bond donors≤5\n0≤H-bond acceptors≤5\n20≤Polar surface area≤130 A˚2'},
                {'subtitle': 'Ghose',                   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'Ghose:\n160≤Molecular weight≤480 daltons\n−0.4≤Log P≤5.6\n40≤Molar refractivity≤130\n20≤Polar surface area≤130 A˚2'},
                {'subtitle': 'GSK Rule',                'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'GSK Rule:\n160≤Molecular weight≤480 daltons\n−0.4≤Log P≤5.6\n0≤H-bond donors≤5\n0≤H-bond acceptors≤5\n20≤Polar surface area≤130 A˚2'},
                {'subtitle': 'Veber (GSK) filter',   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'Veber (GSK) filter: A filter used in drug discovery to evaluate the oral bioavailability of compounds based on the number of rotatable bonds and polar surface area.'},
                {'subtitle': 'Golden Triangle',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'Golden Triangle: A concept in medicinal chemistry describing the balance between lipophilicity, size, and polarity of drug-like compounds.'},
                {'subtitle': 'PAINS',                   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'PAINS: Pan-Assay Interference Compounds, a list of substructures commonly found in screening libraries that may interfere with assay results'},
                {'subtitle': 'Leadlikeness',            'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'Leadlikeness: Criteria used to assess the suitability of a compound as a lead in drug discovery based on its physicochemical properties.'},
                {'subtitle': 'Brenk',                   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': 'Brenk: A set of filters used to identify problematic compounds in drug discovery based on their physicochemical properties.'}
            ]},
            # content : Distribution
            {'title': 'Distribution',
             'contents': [
                {'subtitle': 'PPB',                     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'Fu',                      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'Vd',                      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'BBB Penetration',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'Log Kp',                  'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'Organic anion',           'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'HPP Binding',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'GI absorption ',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'P-gp substrate ',      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
            ]},
            # content : Druglikeness
            {'title': 'Druglikeness',
             'contents': [
                {'subtitle': 'Egan',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'Bioavailability Score',   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
            ]},
            # content : Absorption
            {'title': 'Absorption',
             'contents': [
                {'subtitle': 'Caco-2 Permeability',     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'MDCK Permeability',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'Pgp-inhibitor',           'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'Pgp-substrate',           'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                {'subtitle': 'HIA',                     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
            ]},
            # content : Excretion
            {'title': 'Excretion',
             'contents': [
                 {'subtitle': 'Clearance',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Half-life',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
             ]},
            # content : Metabolism
            {'title': 'Metabolism',
             'contents': [
                 {'subtitle': 'CYP1A2 inhibitor',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'CYP1A2 substrate',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'CYP2C19 inhibitor',      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'CYP2C19 substrate',      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'CYP2C9 inhibitor',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'CYP2C9 substrate',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'CYP2D6 inhibitor',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'CYP2D6 substrate',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'CYP3A4 inhibitor',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'CYP3A4 substrate',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
             ]},
            # content : Toxicity
            {'title': 'Biological Target-Based Toxicity',
             'contents': [
                 {'subtitle': 'hERG Blockers',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'DILI',                       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'AMES Toxicity',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Carcinogencity',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Eye Irritation',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Eye Corrosion',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Hepatotoxicity',             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Rat Oral Acute Toxicity',    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Skin Sensitization',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
             ]},
            # content : Cytotoxicity
            {'title': 'Cytotoxicity',
             'contents': [
                 {'subtitle': 'HepG2',                      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'NIH 3T3',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
             ]},
            # content : Environmental Toxicity
            {'title': 'Environmental Toxicity',
             'contents': [
                 {'subtitle': 'Bioconcentration Factors',   'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Biodegration factors',       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'LD50',                       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
             ]},
            # content : Biological Pathway Toxicity
            {'title': 'Biological Pathway Toxicity',
             'contents': [
                 {'subtitle': 'NR-AR (Androgen receptor)',          'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'NR-AR-LBD (ligand -binding)',        'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'NR-AhR (Aryl hydrocarbon)',          'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'NR-Aromatase',                       'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'NR-ER (Estrogen receptor)',          'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'NR-ER-LBD',                          'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'NR-PPAR-gamma (Peroxisome)',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'SR-ARE (Antioxidant response)',      'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'SR-ATAD5 (ATPase family domain)',    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'SR-HSE (Heat shock factor)',         'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'SR-MMP (Mitochondrial membrane)',    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'SR-p53',                             'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
             ]},
            # content : Toxicophore Rules
            {'title': 'Toxicophore Rules',
             'contents': [
                 {'subtitle': 'Acute Toxicity Rule',                'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Genotoxic Carcinogenicity Rule',     'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'NonGenotoxic Carcinogenicity Rule',  'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Skin Sensitization Rule',            'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'Aquatic Toxicity Rule',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'NonBiodegradable Rule',              'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'SureChEMBL Rule',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '},
                 {'subtitle': 'FAF-Drugs4 Rule',                    'model1': 'num',    'model2': 'num',    'model3': 'num',    'color': '#103B59', 'tooltipText': '   '}
             ]},
        ]

        # 쿼리 파라미터에서 이미지 URL과 SMILES 문자열을 추출
        image_url = request.args.get('image_url')
        smiles = request.args.get('smiles')

        return render_template('evaluation.html', properties=properties, image_url=image_url, smiles=smiles)
    else:
        return 'Molecule not found', 404
