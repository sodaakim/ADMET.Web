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

# endregion

"""
@app.route('/')
def home():
    return 'homeeee'
"""


@app.route('/')
def home():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0',
         'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.',
         'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1',
         'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.',
         'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2',
         'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.',
         'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2',
         'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.',
         'image_path': 'img/path2.jpg'},
    ]
    return render_template('home.html', cards=cards)


@app.route('/kekule')
def Kekule():
    return render_template('box/kekule.html')


@app.route('/example')
def example():
    return render_template('example.html')


if __name__ == '__main__':
    app.run(host='0.0.0.0', port='8000')
# endregion
