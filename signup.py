from flask import Flask, render_template, request, jsonify, redirect, url_for, session
from rdkit import Chem
from datetime import datetime
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
import os
from flask import Blueprint

# Blueprint 객체 생성
signup_route = Blueprint('signup_route', __name__)

@signup_route.route('/signup')
def signupform():
    return render_template('signup.html')