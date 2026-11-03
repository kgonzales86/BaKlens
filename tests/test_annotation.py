import json
from pathlib import Path
from unittest.mock import patch
import pytest
from baklens.annotation import parse_bakta,run_bakta
from baklens.plots import feature_segments

ROOT=Path(__file__).parents[1]

def test_demo_import():
    records,features,meta=parse_bakta((ROOT/'examples/demo.bakta.json').read_bytes())
    assert len(records)==2 and len(features)>0
    assert 'synthetic' in meta['version']['baklens_fixture']

def fixture():
    return {'sequences':[{'id':'a','nt':'ACGTACGT'}],
            'features':[{'sequence':'a','type':'CDS','start':2,'stop':5,'strand':'+'}]}

def test_mismatch():
    with pytest.raises(ValueError,match='match'): parse_bakta(json.dumps(fixture()),{'a':'AAAA'})

def test_invalid_coordinate():
    data=fixture(); data['features'][0]['stop']=9
    with pytest.raises(ValueError,match='outside'): parse_bakta(json.dumps(data))

def test_unknown_sequence():
    data=fixture(); data['features'][0]['sequence']='absent'
    with pytest.raises(ValueError): parse_bakta(json.dumps(data))

def test_origin_spanning_feature():
    data=fixture(); data['features'][0].update(start=7,stop=2,edge=True)
    _,features,_=parse_bakta(json.dumps(data))
    assert [(a,b) for _,a,b in feature_segments(features,8)]==[(7,8),(1,2)]

def test_no_features():
    data=fixture();data['features']=[]
    _,features,_=parse_bakta(json.dumps(data))
    assert features.empty and 'product' in features.columns

@pytest.mark.parametrize('payload',['garbage','[]','{}','{"sequences":[],"features":[]}'])
def test_malformed_json(payload):
    with pytest.raises(ValueError): parse_bakta(payload)

def test_missing_executable():
    with patch('baklens.annotation.shutil.which',return_value=None):
        with pytest.raises(ValueError,match='not installed'): run_bakta({'a':'ACGT'},'/missing')

def test_invalid_db():
    with patch('baklens.annotation.shutil.which',return_value='/bin/bakta'):
        with pytest.raises(ValueError,match='version.json'): run_bakta({'a':'ACGT'},'/missing')

def test_runner_command_and_import(tmp_path):
    db=tmp_path/'db'; db.mkdir(); (db/'version.json').write_text('{}')
    class FakeProcess:
        returncode=0
        def __init__(self,command,**kwargs):
            assert isinstance(command,list) and '--keep-contig-headers' in command
            assert kwargs['start_new_session'] is True
            output=Path(command[command.index('--output')+1]);output.mkdir()
            (output/'baklens.json').write_text(json.dumps(fixture()))
        def wait(self,timeout): return 0
        def __enter__(self): return self
        def __exit__(self,*args): pass
    with patch('baklens.annotation.shutil.which',return_value='/bin/bakta'),patch('baklens.annotation.subprocess.Popen',FakeProcess):
        result,_=run_bakta({'a':'ACGTACGT'},str(db))
        assert json.loads(result)['sequences'][0]['id']=='a'
