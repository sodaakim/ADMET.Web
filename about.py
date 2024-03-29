from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint

# Blueprint 객체 생성
about_route = Blueprint('about_route', __name__)

@about_route.route('/introduce')
def introduce():
    cards = [
        {'title': 'Webserver Advantages',
         'subtitle': 'Sub 0',
         'description': 'Our webserver boasts its capability to predict ADMET profiles quickly and accurately using advanced algorithms and simulation tools as its primary advantage. This feature significantly accelerates the drug development process and reduces unnecessary costs. Additionally, the user-friendly interface and high prediction accuracy ensure that researchers can obtain reliable data on potential drug candidates.',
         'image_path': 'img/path0.jpg'},

        {'title': 'Re-engineered modules and batch evaluation support',
         'subtitle': 'Sub 1',
         'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.',
         'image_path': 'img/path1.jpg'},

        {'title': 'Robust and accurate MGA models',
         'subtitle': 'Sub 2',
         'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.',
         'image_path': 'img/path2.jpg'},

        {'title': 'Practical explanation and guidance',
         'subtitle': 'Sub 2',
         'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.',
         'image_path': 'img/path2.jpg'},
    ]
    return render_template('introduce.html', cards=cards)


@about_route.route('/service_manual')
def service_manual():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('service_manual.html', cards=cards)

@about_route.route('/documentation')
def documentation():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
    ]
    return render_template('documentation.html', cards=cards)


