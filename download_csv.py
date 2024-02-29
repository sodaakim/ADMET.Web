from flask import Flask, Response
import csv
from io import StringIO
from flask import Blueprint

# Blueprint 객체 생성
download_csv_route = Blueprint('download_csv_route', __name__)

@download_csv_route.route('/download-csv', methods=['POST'])
def download_csv():
    '''
    # 세션에서 분석 결과를 가져옵니다.
    if 'results' in session:
        results = session['results']
        # CSV 데이터 생성
        csv_data = [["SMILES", "Formula", "Molecular Weight", "Image URL"]]
        for result in results:
            csv_data.append([result.get('SMILES', 'N/A'), result.get('Formula', 'N/A'),
                             result.get('Molecular Weight', 'N/A'), result.get('Image URL', 'N/A')])

        # CSV 파일로 변환
        si = io.StringIO()
        cw = csv.writer(si)
        cw.writerows(csv_data)
        output = make_response(si.getvalue())
        output.headers["Content-Disposition"] = "attachment; filename=analysis_results.csv"
        output.headers["Content-type"] = "text/csv"
        return output
    else:
        return "No results available", 404
        '''


    # 예제 데이터
    molecule_name = "Example Molecule"
    smiles_notation = "CC(C)OC(=O)C1=CC=CC=C1C(=O)N"
    properties = [
        {"title": "Property 1", "model1": "Value 1", "model2": "Value 2", "model3": "Value 3"},
        {"title": "Property 2", "model1": "Value 4", "model2": "Value 5", "model3": "Value 6"},
        # 여기에 추가적인 속성을 추가할 수 있습니다.
    ]

    # CSV 데이터를 저장할 StringIO 객체 생성
    si = StringIO()
    cw = csv.writer(si)

    # CSV 헤더 작성
    cw.writerow(["Molecule Name", "SMILES Notation", "Property", "Model 1", "Model 2", "Model 3"])

    # 데이터 행 작성
    for prop in properties:
        cw.writerow([molecule_name, smiles_notation, prop["title"], prop["model1"], prop["model2"], prop["model3"]])

    # Flask Response 객체 생성하여 CSV 파일로 전송
    output = si.getvalue()
    return Response(output, mimetype="text/csv",
                    headers={"Content-Disposition": "attachment;filename=molecule_data.csv"})


