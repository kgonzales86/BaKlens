"""Interactive Plotly figures; every plotted coordinate is 1-based."""
import html
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

TEAL = '#0d9488'
BLUE = '#2563eb'
PALETTE = {'cds': TEAL, 'sorf': '#2dd4bf', 'trna': '#f59e0b', 'rrna': '#ec4899',
           'ncrna': '#8b5cf6', 'crispr': '#f97316', 'gap': '#64748b'}


def finish(fig, title, height=400):
    fig.update_layout(title=title, template='plotly_white', height=height,
                      font=dict(family='Arial', color='#172b45'), margin=dict(l=45,r=25,t=65,b=45),
                      paper_bgcolor='white', plot_bgcolor='white')
    return fig


def composition_plot(frame):
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=.08,
                        subplot_titles=('GC content', 'GC skew', 'Ambiguous bases'))
    for row, column, color in [(1,'gc_pct',TEAL),(2,'gc_skew',BLUE),(3,'ambiguous_pct','#f59e0b')]:
        fig.add_trace(go.Scatter(x=frame.position, y=frame[column], mode='lines',
                                name=column, line=dict(color=color), connectgaps=False), row=row,col=1)
    fig.add_hline(y=0, line_dash='dot', line_color='#94a3b8', row=2,col=1)
    fig.update_yaxes(title_text='%', range=[0,100], row=1,col=1)
    fig.update_yaxes(title_text='(G−C)/(G+C)', range=[-1.05,1.05],row=2,col=1)
    fig.update_yaxes(title_text='%', range=[0,100],row=3,col=1)
    fig.update_xaxes(title_text='Position (bp)',row=3,col=1)
    fig.update_layout(showlegend=False)
    return finish(fig, 'Composition along the sequence', 650)


def contig_scatter(table):
    return finish(px.scatter(table, x='length_bp', y='gc_pct', color='ambiguous_pct',
                              hover_name='sequence', log_x=True, color_continuous_scale='Teal',
                              labels={'length_bp':'Length (bp, log scale)', 'gc_pct':'GC (%)',
                                      'ambiguous_pct':'Ambiguous (%)'}), 'Compare contigs')


def length_plot(table):
    data = table.nlargest(30,'length_bp').sort_values('length_bp')
    return finish(px.bar(data,x='length_bp',y='sequence', orientation='h',
                          color_discrete_sequence=[TEAL], labels={'length_bp':'Length (bp)','sequence':''}),
                  'Longest contigs · up to 30', max(350,len(data)*24))


def feature_segments(features, length):
    for row in features.to_dict('records'):
        spans = [(row['start'],row['end'])] if not row['wraps_origin'] else [(row['start'],length),(1,row['end'])]
        for start, end in spans:
            yield row, start, end


def feature_map(features, length, region=None):
    fig = go.Figure()
    for kind, group in features.groupby('type'):
        xs, ys, labels = [], [], []
        for row, start, end in feature_segments(group, length):
            if region and (end < region[0] or start > region[1]):
                continue
            y = 1 if row['strand'] == '+' else -1 if row['strand'] == '-' else 0
            label = html.escape(f"{row['locus']} · {row['gene']} · {row['product']} | {row['start']}–{row['end']} ({row['strand']})")
            xs += [start, end, None]; ys += [y, y, None]; labels += [label,label,None]
        fig.add_trace(go.Scatter(x=xs,y=ys, mode='lines',line=dict(width=9,color=PALETTE.get(kind.lower(),'#64748b')),
                                name=kind,text=labels,hovertemplate='%{text}<extra></extra>'))
    fig.update_yaxes(tickvals=[-1,0,1],ticktext=['Reverse strand','Unstranded','Forward strand'],range=[-1.8,1.8])
    fig.update_xaxes(title='Position (bp)',range=region or [1,length])
    return finish(fig, 'Annotation tracks', 320)


def circular_plot(frame, length, features=None):
    # Reduce displayed windows only; numeric exports retain full resolution.
    sample = frame.iloc[::max(1,int(np.ceil(len(frame)/3000))) ]
    theta = 360 * (sample.position - 1) / length
    fig = go.Figure()
    for name, values, radius, color, scale in [
        ('GC content',sample.gc_pct,2.0,TEAL,.006),
        ('GC skew',sample.gc_skew,1.25,BLUE,.25)]:
        baseline = 50 if name == 'GC content' else 0
        fig.add_trace(go.Scatterpolar(theta=theta,r=radius+(values-baseline)*scale,
            mode='lines',name=name,line=dict(color=color,width=2),customdata=sample[['position']].values,
            text=values.round(3),hovertemplate='%{customdata[0]:,.0f} bp<br>%{text}<extra>%{fullData.name}</extra>'))
    if features is not None:
        for kind, group in features.groupby('type'):
            angles, radii, labels = [], [], []
            for row,start,end in feature_segments(group,length):
                ts=np.linspace(360*(start-1)/length,360*end/length,max(2,int((end-start+1)/length*360)+1))
                radius=2.7 if row['strand']=='+' else 2.5
                label=html.escape(f"{row['locus']} · {row['product']} | {row['start']}–{row['end']} ({row['strand']})")
                angles.extend([*ts,None]); radii.extend([*[radius]*len(ts),None]); labels.extend([*[label]*len(ts),None])
            fig.add_trace(go.Scatterpolar(theta=angles,r=radii,mode='lines',name=kind,
                line=dict(width=5,color=PALETTE.get(kind.lower(),'#64748b')),text=labels,
                hovertemplate='%{text}<extra></extra>'))
    ticks=list(range(0,360,60))
    fig.update_layout(polar=dict(radialaxis=dict(visible=False,range=[0,3]),
        angularaxis=dict(direction='clockwise',rotation=90,tickvals=ticks,
                         ticktext=[f'{round(t/360*length)+1:,} bp' for t in ticks])))
    return finish(fig,'Circular overview · selected sequence',560)
