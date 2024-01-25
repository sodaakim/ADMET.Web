from flask import Flask, request, jsonify, render_template, url_for
from flask_mail import Mail, Message
from itsdangerous import URLSafeTimedSerializer, SignatureExpired

app = Flask(__name__)

# Flask-Mail 설정
app.config['MAIL_SERVER'] = 'smtp.example.com'
app.config['MAIL_PORT'] = 587
app.config['MAIL_USE_TLS'] = True
app.config['MAIL_USERNAME'] = 'your-email@example.com'
app.config['MAIL_PASSWORD'] = 'your-password'
app.config['SECRET_KEY'] = 'your-secret-key'  # 시크릿 키 설정
mail = Mail(app)

# URLSafeTimedSerializer 인스턴스 생성
s = URLSafeTimedSerializer(app.config['SECRET_KEY'])

@app.route('/')
def home():
    cards = [
        {'title': 'Advanced Model for Accurate ADMET', 'subtitle': 'Sub 0', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path0.jpg'},
        {'title': 'Re-engineered modules and batch evaluation support', 'subtitle': 'Sub 1', 'description': 'Based on the studies conducted in papers, we have developed a model capable of accurately predicting ADMET properties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path1.jpg'},
        {'title': 'Robust and accurate MGA models', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        {'title': 'Practical explanation and guidance', 'subtitle': 'Sub 2', 'description': 'Bproperties. Additionally, this model demonstrates excellent performance in predicting various toxicities, proving its robustness and reliability in the field of pharmacological research.', 'image_path': 'img/path2.jpg'},
        # ... 추가 카드 데이터
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

@app.route('/example')
def example():
    return render_template('example.html')

@app.route('/signupform')
def signupform():
    return  render_template('signup.html')

@app.route('/evaluation')
def evaluation():
    return  render_template('evaluation.html')


@app.route('/signup', methods=['POST'])

def signup():
    if request.method == 'POST':

        email = request.json.get('email')
        password = request.json.get('password')

        # Validate email format
        if not re.match(r"[^@]+@[^@]+\.[^@]+", email):
            return jsonify(message="Invalid email format."), 400

        # Hash the password for secure storage
        # hashed_password = generate_password_hash(password)

        # Implement your signup logic here, e.g., save user to database
        # ...

        # Email confirmation token and link
        token = s.dumps(email, salt='email-confirm')
        confirm_url = url_for('confirm_email', token=token, _external=True)

        # Send confirmation email
        msg = Message("Please confirm your email", sender='your-email@example.com', recipients=[email])
        msg.body = f"Please click on the link to confirm your email: {confirm_url}"
        mail.send(msg)

        return jsonify(message="Please check your email to confirm your registration."), 200
    else:
        return 0 #render_template('signup.html')

@app.route('/confirm/<token>')
def confirm_email(token):
    try:
        email = s.loads(token, salt='email-confirm', max_age=3600)
    except SignatureExpired:
        return jsonify(message="The confirmation link is expired."), 400

    # Implement your token confirmation logic here
    # ...

    return jsonify(message="Your email has been confirmed."), 200

if __name__ == '__main__':
    app.run(debug=True)

