import os
from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
from datetime import datetime
from flask_mail import Mail, Message
from itsdangerous import URLSafeTimedSerializer, SignatureExpired
from flask_cors import CORS

import os

app = Flask(__name__)
app.secret_key = 'very_secret_key'  # 실제 배포 시 안전한 키 사용

CORS(app)

@app.route('/')
def home():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('home.html', cards=cards)

@app.route('/introduce')
def introduce():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('introduce.html', cards=cards)

@app.route('/reference')
def reference():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('reference.html', cards=cards)


@app.route('/service_manual')
def service_manual():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('service_manual.html', cards=cards)

@app.route('/api_manual')
def api_manual():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('api_manual.html', cards=cards)


@app.route('/ssbio_lab')
def ssbio_lab():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('ssbio_lab.html', cards=cards)


@app.route('/documentation')
def documentation():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('documentation.html', cards=cards)



@app.route('/kekule')
def Kekule():
    return render_template('box/kekule.html')

@app.route('/analyze', methods=['GET', 'POST'])
def analyze():
    if request.method == 'POST':
        smiles = request.json['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
            filename = f"molecule_{timestamp}.png"
            image_path = os.path.join('static', 'images', filename)
            img = Draw.MolToImage(mol)
            img.save(image_path)
            image_url = url_for('static', filename=f'images/{filename}')

            # 분자 속성 계산
            mol_formula = rdMolDescriptors.CalcMolFormula(mol)
            mol_weight = Descriptors.MolWt(mol)
            mol_volume = Descriptors.MolMR(mol)  # Volume의 대략적 추정
            density = mol_weight / mol_volume if mol_volume else 0  # Density 계산
            num_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)
            num_aromatic_heavy_atoms = len([atom for atom in mol.GetAtoms() if atom.GetIsAromatic()])
            formal_charge = Chem.GetFormalCharge(mol)

            # 세션에 속성 저장
            session['properties'] = {
                'Formula': mol_formula,
                'Molecular Weight': mol_weight,
                'Volume': mol_volume,
                'Density': density,
                'Num. heavy atoms': num_heavy_atoms,
                'Num. arom. heavy atoms': num_aromatic_heavy_atoms,
                'formal Charge': formal_charge
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

@app.route('/screening')
def screening():
    steps = [
        {'title': 'Step 4', 'description': 'paste admet1'},
        {'title': 'Step 5', 'description': 'paste admet2'},
        {'title': 'Step 6', 'description': 'paste admsdsdsdsdsd'}
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
    return render_template('screening.html', steps=steps, cardsA=cardsA, cardsI=cardsI)

@app.route('/example')
def example():
    return render_template('example.html')

@app.route('/signup')
def signupform():
    return render_template('signup.html')

@app.route('/evaluation')
def evaluation():
    smiles = request.args.get('smiles')
    if 'properties' in session:
        prop_values = session['properties']  # analyze 함수에서 계산 및 저장한 값
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

            # content : Biological Target-Based Toxicity
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
                 {'subtitle': 'Eye Irritation',         'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Eye Irritation',         'num': 'num', 'color': '#103B59'},
                 {'subtitle': 'Eye Irritation',         'num': 'num', 'color': '#103B59'},
             ]},
            # content : Functional Toxicity
            {'title': 'Functional Toxicity',
             'contents': [
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

@app.route('/result')
def result():
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
            {'subtitle': 'Formula',                 'num': 'num', 'color': '#103B59'},
            {'subtitle': 'Molecular Weight',        'num': 'num', 'color': '#103B59'},
            {'subtitle': 'Volume',                  'num': 'num', 'color': '#103B59'},
            {'subtitle': 'Density',                 'num': 'num', 'color': '#103B59'},
            {'subtitle': 'Num. heavy atoms',        'num': 'num', 'color': '#103B59'},
            {'subtitle': 'Num. arom. heavy atoms',  'num': 'num', 'color': '#103B59'},
            {'subtitle': 'formal Charge',           'num': 'num', 'color': '#103B59'},
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
        {'title': 'Medicinal Chemistry', 'color': cardColors[1],
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
        {'title': 'Distribution', 'color': cardColors[2],
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
        {'title': 'Druglikeness', 'color': cardColors[3],
         'contents': [
            {'subtitle': 'Egan',                    'num': 'num', 'color': '#103B59'},
            {'subtitle': 'Bioavailability Score',   'num': 'num', 'color': '#103B59'}
        ]},
        # content : Absorption
        {'title': 'Absorption', 'color': cardColors[4],
         'contents': [
            {'subtitle': 'Caco-2 Permeability',     'num': 'num', 'color': '#103B59'},
            {'subtitle': 'MDCK Permeability',       'num': 'num', 'color': '#103B59'},
            {'subtitle': 'Pgp-inhibitor',           'num': 'num', 'color': '#103B59'},
            {'subtitle': 'Pgp-substrate',           'num': 'num', 'color': '#103B59'},
            {'subtitle': 'HIA',                     'num': 'num', 'color': '#103B59'}
        ]},
        # content : Excretion
        {'title': 'Excretion', 'color': cardColors[5],
         'contents': [
             {'subtitle': 'Clearance',              'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Half-life',              'num': 'num', 'color': '#103B59'}
         ]},
        # content : Metabolism
        {'title': 'Metabolism', 'color': cardColors[6],
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

        # content : Biological Target-Based Toxicity
        {'title': 'Biological Target-Based Toxicity', 'color': cardColors[7],
         'contents': [
             {'subtitle': 'hERG Blockers',          'num': 'num', 'color': '#103B59'},
             {'subtitle': 'DILI',                   'num': 'num', 'color': '#103B59'},
             {'subtitle': 'AMES Toxicity',          'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Carcinogencity',         'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Eye Irritation',         'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Eye Corrosion',          'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Hepatotoxicity',         'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Rat Oral Acute Toxicity', 'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Eye Irritation',         'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Eye Irritation',         'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Eye Irritation',         'num': 'num', 'color': '#103B59'},
         ]},
        # content : Functional Toxicity
        {'title': 'Functional Toxicity', 'color': cardColors[8],
         'contents': [
             {'subtitle': 'Skin Sensitization',     'num': 'num', 'color': '#103B59'}
         ]},
        # content : Cytotoxicity
        {'title': 'Cytotoxicity', 'color': cardColors[9],
         'contents': [
             {'subtitle': 'HepG2',                  'num': 'num', 'color': '#103B59'},
             {'subtitle': 'NIH 3T3',                'num': 'num', 'color': '#103B59'}
         ]},
        # content : Environmental Toxicity
        {'title': 'Environmental Toxicity', 'color': cardColors[10],
         'contents': [
             {'subtitle': 'Bioconcentration Factors',   'num': 'num', 'color': '#103B59'},
             {'subtitle': 'Biodegration factors',       'num': 'num', 'color': '#103B59'},
             {'subtitle': 'LD50',                       'num': 'num', 'color': '#103B59'}
         ]},
        # content : Biological Pathway Toxicity
        {'title': 'Biological Pathway Toxicity', 'color': cardColors[11],
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
        {'title': 'Toxicophore Rules', 'color': cardColors[12],
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

    return render_template('result.html', properties=properties)



if __name__ == '__main__':
    app.run(debug=True, port=8000)



