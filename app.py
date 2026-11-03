"""Run: python -m streamlit run app.py"""
from pathlib import Path
from io import BytesIO
import hashlib
import json
import os
import zipfile
import streamlit as st
import plotly.express as px
from baklens.analysis import parse_fasta, sequence_table, assembly_stats, windows, summary_markdown
from baklens.annotation import parse_bakta, run_bakta, bakta_available
from baklens.plots import composition_plot, contig_scatter, length_plot, feature_map, circular_plot, finish

st.set_page_config(page_title='BaKlens | Genome explorer',page_icon='🧬',layout='wide')
st.markdown('''<style>
.block-container {padding-top:2.4rem; max-width:1450px;}
[data-testid="stMetric"] {background:white;border:1px solid #e2e8f0;border-radius:14px;padding:18px;}
h1 {letter-spacing:-1.8px;} .eyebrow {color:#0d9488;font-size:12px;letter-spacing:3px;font-weight:700;}
</style>''',unsafe_allow_html=True)
st.markdown('<div class="eyebrow">BACTERIAL GENOMICS / LOCAL EXPLORATION</div>',unsafe_allow_html=True)
st.title('BaKlens')
st.write('From DNA sequence to a clearer view of your genome.')
st.caption('Explore composition, inspect contigs, and connect Bakta annotations. No API key required.')

with st.sidebar:
    st.header('Your workspace')
    mode=st.radio('Data source',['Synthetic demo','DNA FASTA','Bakta JSON'])
    upload=None
    if mode=='DNA FASTA':
        upload=st.file_uploader('Upload bacterial DNA',type=['fasta','fa','fna','txt'])
    elif mode=='Bakta JSON':
        upload=st.file_uploader('Upload Bakta results',type=['json'])
    st.divider()
    st.caption('Files are processed by the machine running this app. The local launch stays on your computer. No external uploads or AI calls are made.')
    st.caption('Version 0.1.0 · Research and education')

features=None
provenance={}
raw_json=None
try:
    if mode=='Synthetic demo':
        raw_json=(Path(__file__).parent/'examples/demo.bakta.json').read_bytes()
        records,features,provenance=parse_bakta(raw_json)
        st.info('Synthetic demonstration: these sequences and illustrative features are generated fixtures, not real Bakta predictions.')
    elif upload is None:
        st.info('Choose a FASTA or Bakta JSON file in the sidebar to start. Try the synthetic demo to explore the interface.')
        st.stop()
    elif mode=='DNA FASTA':
        records=parse_fasta(upload.getvalue().decode('utf-8-sig'))
    else:
        raw_json=upload.getvalue()
        records,features,provenance=parse_bakta(raw_json)
except (ValueError,UnicodeDecodeError) as exc:
    st.error(str(exc)); st.stop()

fingerprint=hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest()
if mode=='DNA FASTA':
    with st.expander('Add Bakta annotations',expanded=False):
        st.write('Import results, or run Bakta locally after installing its software and database. Bacterial genomes and plasmids only.')
        imported=st.file_uploader('Matching Bakta JSON',type=['json'],key='annotations')
        if imported:
            try:
                raw_json=imported.getvalue()
                _,features,provenance=parse_bakta(raw_json,records)
            except ValueError as exc:
                raw_json=None; st.error(str(exc))
        st.caption('Local annotation can take minutes or longer. Uploads are not sent to Bakta Web.')
        if not bakta_available():
            st.info('Bakta executable not found. FASTA plots and JSON import work without it. See README for installation.')
        db=os.environ.get('BAKTA_DB','')
        st.caption('The local database is configured using the BAKTA_DB environment variable.')
        threads=st.number_input('Annotation threads',min_value=1,max_value=32,value=2)
        if st.button('Run local Bakta',disabled=not bakta_available() or not db):
            try:
                with st.spinner('Bakta is annotating the genome. This may take several minutes…'):
                    result,logs=run_bakta(records,db,int(threads))
                st.session_state['annotation_result']=(fingerprint,result,logs)
            except ValueError as exc:
                st.error(str(exc))
        saved=st.session_state.get('annotation_result')
        if saved and saved[0]==fingerprint and not imported:
            raw_json=saved[1]
            _,features,provenance=parse_bakta(raw_json,records)
            st.success('Local Bakta annotation complete.')
            st.download_button('Download annotation log',saved[2],file_name='baklens-bakta.log')

stats=assembly_stats(records)
cols=st.columns(4)
cols[0].metric('TOTAL BASES',f"{stats['total_bp']:,}")
cols[1].metric('SEQUENCES',f"{stats['sequences']:,}")
cols[2].metric('GC CONTENT',f"{stats['gc_pct']:.2f}%")
cols[3].metric('N50',f"{stats['n50_bp']:,} bp")
st.caption(f"L50: {stats['l50']:,} · Ambiguous bases: {stats['ambiguous_bp']:,} · GC% uses A/C/G/T bases only.")

with st.sidebar:
    st.divider()
    selected=st.selectbox('Sequence',list(records))
    size=int(st.number_input('Window size (bp)',min_value=1,max_value=1_000_000,value=1000,step=100))
    step=int(st.number_input('Step (bp)',min_value=1,max_value=1_000_000,value=1000,step=100))
    st.caption('Windows stop at sequence boundaries. A final partial window is retained.')
seq=records[selected]
try:
    frame=windows(seq,size,step)
except ValueError as exc:
    st.warning(str(exc)); st.stop()
table=sequence_table(records)
selected_features=features[features.sequence==selected] if features is not None else None
report=summary_markdown(records,size,step,features,json.dumps(provenance,indent=2))
fig_comp=composition_plot(frame)
fig_circle=circular_plot(frame,len(seq),selected_features)
fig_scatter=contig_scatter(table)
fig_lengths=length_plot(table)

overview,composition,genome,annotation,exports=st.tabs(['Assembly overview','GC & skew','Genome map','Annotations','Export & methods'])
with overview:
    left,right=st.columns([1.1,1])
    with left: st.plotly_chart(fig_scatter,use_container_width=True)
    with right: st.plotly_chart(fig_lengths,use_container_width=True)
    st.dataframe(table,use_container_width=True,hide_index=True)
    st.caption('Unusual composition is a prompt for further investigation, not a contamination diagnosis. N50 measures contiguity, not correctness.')
with composition:
    st.subheader(f'Composition · {selected}')
    st.plotly_chart(fig_comp,use_container_width=True)
    st.caption('GC skew = (G − C) / (G + C). Positive values mean more G than C on the supplied strand. Undefined values appear as gaps. Skew alone does not identify a replication origin.')
with genome:
    st.plotly_chart(fig_circle,use_container_width=True)
    st.caption('Layout only: circular plotting does not imply biological circularity. Rings: outer forward features, reverse/unstranded features, GC% centred at 50%, and GC skew centred at zero. Display is sampled to at most ~3,000 windows.')
    if selected_features is not None and not selected_features.empty:
        left,right=st.columns(2)
        start=int(left.number_input('Region start (bp)',1,len(seq),1,key=f'start-{fingerprint}-{selected}'))
        end=int(right.number_input('Region end (bp)',1,len(seq),len(seq),key=f'end-{fingerprint}-{selected}'))
        if start<=end:
            st.plotly_chart(feature_map(selected_features,len(seq),[start,end]),use_container_width=True)
        else: st.warning('Region start must not exceed region end.')
    else: st.info('Add Bakta annotations to display gene and RNA tracks.')
with annotation:
    if features is None:
        st.info('No annotations loaded. Use “Add Bakta annotations” above or open a Bakta JSON file.')
    else:
        if not features.empty:
            counts=features.type.value_counts().rename_axis('type').reset_index(name='count')
            st.plotly_chart(finish(px.bar(counts,x='type',y='count',color_discrete_sequence=['#0d9488']),
                                  'Annotated feature counts'),use_container_width=True)
            cds=features[features.type.str.lower()=='cds']
            hypothetical=cds['product'].str.contains('hypothetical protein',case=False,regex=False).sum()
            st.caption(f'{len(cds):,} CDS · {hypothetical:,} labelled hypothetical protein. Other labels are not experimental validation.')
        query=st.text_input('Search gene, locus, product, or sequence')
        filtered=features
        if query:
            mask=features[['gene','locus','product','sequence']].astype(str).apply(lambda c:c.str.contains(query,case=False,regex=False)).any(axis=1)
            filtered=features[mask]
        st.dataframe(filtered,use_container_width=True,hide_index=True)
        with st.expander('Annotation provenance'): st.json(provenance)
with exports:
    st.markdown(report)
    st.download_button('Download summary',report,file_name='baklens-report.md')
    st.download_button('Download contig statistics',table.to_csv(index=False),file_name='baklens-contigs.csv')
    st.download_button('Download selected sequence windows',frame.to_csv(index=False),file_name='baklens-windows.csv')
    if features is not None:
        st.download_button('Download features',features.to_csv(index=False),file_name='baklens-features.csv')
    if raw_json:
        st.download_button('Download Bakta JSON',raw_json,file_name='baklens-annotations.json')
    if st.button('Prepare interactive plot bundle'):
        buffer=BytesIO()
        with zipfile.ZipFile(buffer,'w',zipfile.ZIP_DEFLATED) as archive:
            for name,fig in [('composition',fig_comp),('circular',fig_circle),('contigs',fig_scatter),('lengths',fig_lengths)]:
                archive.writestr(f'{name}.html',fig.to_html(include_plotlyjs=True,full_html=True))
            archive.writestr('report.md',report)
            archive.writestr('contigs.csv',table.to_csv(index=False))
            archive.writestr('windows.csv',frame.to_csv(index=False))
        st.download_button('Download plot bundle',buffer.getvalue(),file_name='baklens-plots.zip')
    st.caption('HTML plots include Plotly and work offline. Use the camera icon on individual plots to save PNGs. Selected-sequence windows are not whole-assembly exports.')
st.divider()
st.caption('BaKlens · Built with Python, Biopython, pandas, Plotly and Streamlit. Optional annotations by Bakta. No generative AI runs inside this application.')
