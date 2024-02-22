from flask import Flask, render_template, request, make_response
from weasyprint import HTML
import json
from flask import Blueprint

# Blueprint 객체 생성
download_pdf_route = Blueprint('download_pdf_route', __name__)

@download_pdf_route.route('/download-pdf', methods=['POST'])
def download_pdf():
    # AJAX 요청으로부터 데이터를 받아옴
    data = request.json
    # pdf_template.html 템플릿을 렌더링하며 데이터를 전달
    rendered_html = render_template('pdf_template.html', **data)
    # 렌더링된 HTML을 PDF로 변환
    pdf = HTML(string=rendered_html).write_pdf()
    # PDF 파일로 응답 생성
    response = make_response(pdf)
    response.headers['Content-Type'] = 'application/pdf'
    response.headers['Content-Disposition'] = 'attachment; filename="report.pdf"'
    return response

