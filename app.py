from flask import Flask, request, jsonify, render_template, url_for
from flask_mail import Mail, Message
from itsdangerous import URLSafeTimedSerializer, SignatureExpired

app = Flask(__name__)
@app.route('/')
def home():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('home.html', cards=cards)

@app.route('/Kekule')
def Kekule():
    return render_template('box/kekule.html')

@app.route('/analyze')
def analyze():
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

@app.route('/signupform')
def signupform():
    return render_template('signup.html')

@app.route('/evaluation')
def evaluation():
    properties = [ # property
        # content : Physicochemical Property
        {'title': 'Physicochemical Properties',
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
        # content : Distribution
        {'title': 'Title6',
         'contents': [
            {'subtitle': 'AA', 'num': 'AAA', 'color': '#103B59'},
            {'subtitle': 'BB', 'num': 'BBB', 'color': '#103B59'},
            {'subtitle': 'BB', 'num': 'BBB', 'color': '#103B59'}
        ]}
    ]
    return render_template('evaluation.html', properties=properties)


if __name__ == '__main__':
    app.run(debug=True)

