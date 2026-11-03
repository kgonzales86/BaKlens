from pathlib import Path
from streamlit.testing.v1 import AppTest

ROOT=Path(__file__).parents[1]

def test_demo_and_navigation():
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=30).run()
    assert not app.exception
    assert app.title[0].value=='BaKlens'
    assert app.metric[0].value=='104,000'
    app.sidebar.selectbox[0].select('demo_contig_2').run()
    assert not app.exception
    app.sidebar.number_input[0].set_value(500).run()
    assert not app.exception
    app.text_input[0].set_value('hypothetical').run()
    assert not app.exception

def test_empty_upload_state():
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=30).run()
    app.sidebar.radio[0].set_value('DNA FASTA').run()
    assert not app.exception
    assert any('Choose a FASTA' in x.value for x in app.info)
