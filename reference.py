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
    data_models = [
        {'title': 'Model Example1', 'image_url': 'static/resources/reference.png', 'description': 'This model is about..', 'link': 'Link'},
        {'title': 'Model Example2', 'image_url': 'static/resources/reference.png', 'description': 'This model is about..', 'link': 'Link'}
    ]

    parameters = [
        {'title': 'parameters Example1', 'image_url': 'static/resources/reference.png', 'description': 'parameters..', 'link': 'Link'},
        {'title': 'parameters Example2', 'image_url': 'static/resources/reference.png', 'description': 'parameters..', 'link': 'Link'}
    ]
    return render_template('reference.html', data_models=data_models, parameters=parameters)


@reference_route.route('/ssbio_lab')
def ssbio_lab():
    lab_projects = [
        {'title': 'lab projects Example1', 'image_url': 'static/resources/reference.png', 'description': 'This model is about..', 'link': 'Link'},
        {'title': 'lab projects Example2', 'image_url': 'static/resources/reference.png', 'description': 'This model is about..', 'link': 'Link'}
    ]
    research_themes = [
        {'title': 'research themes Example1', 'image_url': 'static/resources/reference.png', 'description': 'This model is about..', 'link': 'Link'},
        {'title': 'research themes Example2', 'image_url': 'static/resources/reference.png', 'description': 'This model is about..', 'link': 'Link'}
    ]
    return render_template('ssbio_lab.html', lab_projects=lab_projects, research_themes=research_themes)

