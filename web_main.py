# -*- coding: utf-8 -*-
import os
from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from flask_mail import Mail, Message
from itsdangerous import URLSafeTimedSerializer, SignatureExpired

# region main
app = Flask(__name__)
app.secret_key = 'secret_key'
app.jinja_env.filters['zip'] = zip

# region import other Python Code
from analyze import analyze_route
app.register_blueprint(analyze_route)

from evaluation import evaluation_route
app.register_blueprint(evaluation_route)

from result import result_route
app.register_blueprint(result_route)

from about import about_route
app.register_blueprint(about_route)

from reference import reference_route
app.register_blueprint(reference_route)

from api import api_route
app.register_blueprint(api_route)

from signup import signup_route
app.register_blueprint(signup_route)

from screening import screening_route
app.register_blueprint(screening_route)

from download_pdf import download_pdf_route
app.register_blueprint(download_pdf_route)

from evaluation_result import evaluation_result_route
app.register_blueprint(evaluation_result_route)

# endregion

"""
@app.route('/')
def home():
    return 'homeeee'
"""


@app.route('/')
def home():
    cards = [
        {'title': 'Webserver Advantages',
         'subtitle': 'Sub 0',
         'icon': 'fas fa-poll',
         'description': 'Our webserver boasts its capability to predict ADMET profiles quickly and accurately using advanced algorithms and simulation tools as its primary advantage. This feature significantly accelerates the drug development process and reduces unnecessary costs. Additionally, the user-friendly interface and high prediction accuracy ensure that researchers can obtain reliable data on potential drug candidates.',
         'image_path': 'static/resources/service_manual.png'},

        {'title': 'Expanded ADMET Profiles',
         'subtitle': 'Sub 1',
         'icon': 'fas fa-poll',
         'description': 'The webserver offers an extended ADMET profile covering 88 related characteristics across 7 different categories. This surpasses competitive sites in terms of the number of properties and improved performance. Users can leverage this comprehensive dataset to set their own criteria for promising and desirable molecules.',
         'image_path': 'static/resources/service_manual.png'},

        {'title': 'Enhanced User Experience',
         'subtitle': 'Sub 2',
         'icon': 'fas fa-poll',
         'description': 'Functional modules have been redesigned and optimized to enhance the user experience. Notably, the server supports batch uploading and downloading, allowing users to adjust download options to their preference.',
         'image_path': 'static/resources/screening.png'},

        {'title': 'Reinforced Prediction Models',
         'subtitle': 'Sub 3',
         'icon': 'fas fa-poll',
         'description': 'The simultaneous development of classification and regression predictors, along with the application of deep learning, makes multitask learning seamless, leading to performance enhancement across many modeled endpoints. This advancement enables users to gain an accurate overall ADMET profile for input molecules.',
         'image_path': 'static/resources/service_manual.png'},

        {'title': 'Intuitive Usage',
         'subtitle': 'Sub 4',
         'icon': 'fas fa-poll',
         'description': 'The usage of our webserver is remarkably intuitive. Users simply need to input the SMILES string and click the \'Start Analysis\' button to receive detailed predictive results for the drug\'s ADMET profile. These results include the drug\'s absorption rate, tissue distribution, metabolic pathways, excretion methods, and potential toxicity effects, ensuring a comprehensive understanding of the candidate molecule\'s properties.',
         'image_path': 'static/resources/screening.png'},
    ]
    return render_template('home.html', cards=cards)


@app.route('/kekule')
def Kekule():
    return render_template('box/kekule.html')

@app.route('/pdf')
def PDF():
    return render_template('pdf_template.html')


@app.route('/example')
def example():
    return render_template('example.html')


if __name__ == '__main__':
    app.run(host='0.0.0.0', port='8000')
# endregion
