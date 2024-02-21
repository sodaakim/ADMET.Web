from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint

# Blueprint 객체 생성
reference_route = Blueprint('reference_route', __name__)
@reference_route.route('/reference')
def reference():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('reference.html', cards=cards)


@reference_route.route('/ssbio_lab')
def ssbio_lab():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('ssbio_lab.html', cards=cards)

