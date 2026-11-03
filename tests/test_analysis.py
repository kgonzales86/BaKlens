import math
import pytest
from baklens.analysis import parse_fasta,composition,assembly_stats,windows,summary_markdown

@pytest.mark.parametrize('text',['ACGT','>a\n','>a\nACGT\n>a\nAA','>a\nACGU','>a\nACGT-','>a\nMPEPTIDE'])
def test_bad_fasta(text):
    with pytest.raises(ValueError): parse_fasta(text)

def test_lowercase_iupac():
    assert parse_fasta('>a description\nacgtryn\n')=={'a':'ACGTRYN'}

def test_composition_denominators():
    result=composition('GGGCATNNR')
    assert result['gc_pct']==pytest.approx(100*4/6)
    assert result['gc_skew']==.5
    assert result['ambiguous_bp']==3
    assert result['n_bp']==2

def test_undefined_is_missing():
    assert math.isnan(composition('NNNN')['gc_pct'])
    assert math.isnan(composition('AAAA')['gc_skew'])

def test_weighted_gc_and_n50():
    stats=assembly_stats({'a':'G'*6,'b':'A'*3,'c':'NN'})
    assert stats['n50_bp']==6
    assert stats['l50']==1
    assert stats['gc_pct']==pytest.approx(100*6/9)

def test_partial_window_and_coordinates():
    frame=windows('GGCCAAATTT',4,4)
    assert frame.start.tolist()==[1,5,9]
    assert frame.end.tolist()==[4,8,10]
    assert frame.length_bp.tolist()==[4,4,2]

def test_overlap_stops_when_end_reached():
    frame=windows('ACGTACGTAA',6,3)
    assert frame.start.tolist()==[1,4,7]
    assert frame.end.tolist()==[6,9,10]

def test_short_sequence():
    assert len(windows('ACG',1000,1000))==1

@pytest.mark.parametrize('size,step',[(0,1),(1,0),(-1,1)])
def test_invalid_windows(size,step):
    with pytest.raises(ValueError): windows('ACGT',size,step)

def test_limit_windows():
    with pytest.raises(ValueError): windows('A'*50001,1,1)

def test_report_all_ambiguous():
    report=summary_markdown({'x':'NNN'},10,10)
    assert 'undefined (no A/C/G/T)' in report
    assert 'not generative AI' in report
