import pandas as pd
import numpy as np
import os
from snp_analysis_tools_sherlock import *
from coalescence_analysis_tools import *
import iqplot
import bokeh.plotting
import bokeh.io
import holoviews as hv
from holoviews import dim, opts
import bokeh.models
from bokeh.layouts import gridplot
from scipy.optimize import curve_fit 
pd.options.mode.chained_assignment = None 
hv.extension('bokeh')

tick_font_size ='22px'
label_font_size = '28px'

def get_pred_sel(ts,sel_coeff, start_value,):
    pred=np.exp(ts*sel_coeff)
    pred = pred*start_value/(1 - start_value + start_value*pred)
    return pred 

def export_plot_pdf(p, fname_prefix):
    p.output_backend = "svg"
    plot_fname_svg=f'plots/{fname_prefix}.svg'
    plot_fname_pdf=f'plots/{fname_prefix}.pdf'
    bokeh.io.export_svg(p,filename=plot_fname_svg)
    convert_scripts = f'rsvg-convert -f pdf -o {plot_fname_pdf} {plot_fname_svg}' 
    os.system(convert_scripts)

#df_metadata.head()

def plot_abundance_trajectories_with_in_abundance(freq_df, df_abundance, colorby='meso',alpha=1.,strain1=True,
                                                  shift=True,markers='meso',
                             p1=None,p2=None,emph=[]):
    if colorby=='media':
        cmap = {'fecal': bokeh.palettes.Set2[8][-2],
        'mBHI':bokeh.palettes.Set2[8][1],
        'mGAM': bokeh.palettes.Set2[8][0],}
  #  marker_dic =['circle_cross','circle_dot','circle_x','circle_y','circle']
    marker_dic = ['circle','triangle','square', 'diamond', 'hexagon']

    if not p2:
        
        p2 = bokeh.plotting.figure(width = 500, height = 300)
        p2.xaxis.axis_label_text_font_size=label_font_size
        p2.xaxis.major_label_text_font_size=tick_font_size
        p2.yaxis.axis_label_text_font_size=label_font_size
        p2.yaxis.major_label_text_font_size='10px'
        p2.yaxis.minor_tick_line_color= None
    if not p1:
        p1 = bokeh.plotting.figure(width = 500, height = 150)
        p1.xaxis.axis_label_text_font_size=label_font_size
        p1.xaxis.major_label_text_font_size=tick_font_size
        p1.yaxis.axis_label_text_font_size=label_font_size
        p1.yaxis.major_label_text_font_size=tick_font_size
        p1.yaxis.minor_tick_line_color= None

    for i, type_meso in enumerate(freq_df['type_mesocosm'].unique()):
        mesos = freq_df.loc[freq_df['type_mesocosm'] ==type_meso, 'mesocosm'].unique()
        media = type_meso.split('-')[-1]
        if colorby=='media':
            color = cmap[media]
        if colorby=='meso':
            if media == 'mBHI':
                #cmap = [bokeh.palettes.Light[9][1], bokeh.palettes.Light[9][2], bokeh.palettes.Light[9][3], bokeh.palettes.Bright[7][2]]
               cmap = bokeh.palettes.Light[9][3:]
            else:
                cmap= [bokeh.palettes.Light[9][0], bokeh.palettes.Light[9][5], bokeh.palettes.Light[9][6], bokeh.palettes.Light[9][7]]
                cmap = bokeh.palettes.Light[9][3:]

        for j,mesocosm in enumerate(mesos):
            df_meso = freq_df.loc[freq_df['mesocosm'] == mesocosm,:]
            df_abundance_meso = df_abundance.loc[df_abundance['mesocosm']==mesocosm,:]
            in_sample = df_meso['inoculumn_sample'].unique()[0]
            rep = df_meso['replicate'].unique()[0]
            marker = 'circle'
            if colorby=='meso':
                color=cmap[int(rep)]
            if markers == 'meso':
                marker = marker_dic[int(rep)]

            if in_sample not in df_meso['sample'].values:
                df_meso= pd.concat([df_meso,freq_df.loc[freq_df['sample'] == in_sample,:]],axis=0)
                df_abundance_meso = pd.concat([df_abundance_meso, df_abundance.loc[df_abundance['sample']==in_sample]],axis=0)
            df_meso = df_meso.sort_values(by = 'passage') 
            df_abundance_meso = df_abundance_meso.sort_values(by='passage')
            if strain1:
                minor_strain = 'boot_med1'
                upper = 'boot_high1'
                lower = 'boot_low1'
                if shift:
                    df_meso['boot_med1_shift']=df_meso['boot_med1']/(df_meso['boot_med1']+df_meso['boot_med2'])
                    df_meso['boot_high1_shift']=df_meso['boot_high1']/(df_meso['boot_high1']+df_meso['boot_high2'])
                    df_meso['boot_low1_shift']=df_meso['boot_low1']/(df_meso['boot_low1']+df_meso['boot_low2'])
                    minor_strain = 'boot_med1_shift'
                    upper = 'boot_high1_shift'
                    lower = 'boot_low1_shift'
            else:
                minor_strain = 'boot_med2'
                upper = 'boot_high2'
                lower = 'boot_low2'
            df_source = bokeh.models.ColumnDataSource(data=dict(base=df_meso['passage'].values, 
                                                                upper=df_meso[upper].values, lower=df_meso[lower].values))
            error = bokeh.models.Whisker(base='base', upper='upper', lower='lower', source=df_source,line_color='black',
                                         line_width=1,  level="annotation",)
            error.upper_head.size=5
            error.lower_head.size=5
            p1.line(df_abundance_meso['passage'].values,
                   df_abundance_meso['relative_abundance'].values,color=color,alpha=alpha)
            p1.scatter(df_abundance_meso['passage'].values,
                   df_abundance_meso['relative_abundance'].values,color = color, legend_label = mesocosm, size = 8,
                   marker=marker,#line_color='black',
                      alpha=alpha)
            
            p2.line(df_meso['passage'].values,
                   df_meso[minor_strain].values,color=color,alpha=alpha)
            p2.scatter(df_meso['passage'].values,
                   df_meso[minor_strain].values,color = color, legend_label = mesocosm, size = 8,marker=marker,
                   #line_color='black',
                      alpha=alpha)
            p2.add_layout(error,)
           # bokeh.io.show(p2)
        
    subs_dic = {'AA':'A','AE':'B','AF':'C'}
    sub1, sub2 = type_meso.split('-')[0],type_meso.split('-')[1]
    sub1 = subs_dic[sub1] 
    sub2 = subs_dic[sub2] 
    if strain1:
        p2.yaxis.axis_label = f'Freq Strain {sub1} vs {sub2}'

    else:
        p2.yaxis.axis_label = f'Freq Strain {sub2} vs {sub1}'

    p1.xaxis.axis_label = 'Timepoint'
    # p2.y_range = bokeh.models.Range1d(-.005,1.005)
    # p2.x_range = bokeh.models.Range1d(-.5,7.5)
    p1.xgrid.grid_line_color = None
    p1.ygrid.grid_line_color = None

    p1.yaxis.axis_label = 'Sp Abun'
    p2.xaxis.axis_label = 'Timepoint'
    # p2.y_range = bokeh.models.Range1d(-.005,1.005)
    # p2.x_range = bokeh.models.Range1d(-.5,7.5)
    p2.xgrid.grid_line_color = None
    p2.ygrid.grid_line_color = None
    p2.yaxis.minor_tick_line_color= None


    return p1,p2

def make_plain_linear_plot_abundance(df_abundance, species, inoculumn,emph=[],same_plots = False,colorby='media',markers='meso',):
    e003_metadata = pd.read_csv('e003_coalescence_metadata_round4_good.csv').set_index('sample')
   # print(e003_metadata.columns.values)

    fname1=f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/{inoculumn}_parent1_info.csv'
    df_both = get_both_dfs(fname1)
    df_both['boot_med1_shift'] = df_both['boot_med1']/(df_both['boot_med1']+df_both['boot_med2'])
    df_both = df_both.loc[np.intersect1d(df_both.index.values,e003_metadata.index.values),:]
    df_both_meta = pd.concat([df_both,e003_metadata.loc[np.intersect1d(df_both.index.values,
                                                                               e003_metadata.index.values),:]],axis=1).reset_index()
    df_both_meta_good = df_both_meta.loc[df_both_meta['total_shift']<.1,:]
    df_abundance_good = df_abundance.loc[df_abundance['species_id']==species,:]
    plots_abundance = []
    plots_strain=[]
    type_mesos=[]
    if same_plots:
        p1 = bokeh.plotting.figure(width = 500, height = 190)
        p1.xaxis.axis_label_text_font_size=label_font_size
        p1.xaxis.major_label_text_font_size=tick_font_size
        p1.yaxis.axis_label_text_font_size=label_font_size
        p1.yaxis.major_label_text_font_size=tick_font_size
        p1.yaxis.minor_tick_line_color= None
        p2 = bokeh.plotting.figure(width = 500, height = 300)
        p2.xaxis.axis_label_text_font_size=label_font_size
        p2.xaxis.major_label_text_font_size=tick_font_size
        p2.yaxis.axis_label_text_font_size=label_font_size
        #p2.yaxis.major_label_text_font_size='10px'
        p2.yaxis.minor_tick_line_color= None
        
    current_max = 0
    for good_meso in df_both_meta_good['type_mesocosm'].unique():
        df_both_meta_good_meso= df_both_meta_good.loc[(df_both_meta_good['type_mesocosm'] == good_meso)+(df_both_meta_good['is_inoculumn']),:].copy()
        df_abundance_good_meso = df_abundance_good.loc[(df_abundance_good['type_mesocosm'] == good_meso)+(df_abundance_good['is_inoculumn']),:].copy()
        df_abundance_good_meso['type_mesocosm'] = good_meso
        df_both_meta_good_meso['type_mesocosm'] = good_meso
        if same_plots:
            p1,p2 = plot_abundance_trajectories_with_in_abundance(df_both_meta_good_meso,df_abundance_good_meso, strain1=True,alpha=1,
                                     emph=emph,p1=p1,p2=p2,colorby=colorby,markers=markers)
        else:
            p1,p2 = plot_abundance_trajectories_with_in_abundance(df_both_meta_good_meso, df_abundance_good_meso, strain1=True,alpha=1,
                                     emph=emph,colorby=colorby,markers=markers)
    
        species = fname1.split('/')[-2]
       # species_name = df_metadata.loc[species,'species_plot']
        subject1 = fname1.split('/')[-1].split('-')[0]
        in_name = fname1.split('/')[-1].split('_')[0]
       # p2.title.text = f'{species_name}, {subject1} in {in_name} collisions'
        p2.x_range = bokeh.models.Range1d(-.1,7.1)
        p2.y_range = bokeh.models.Range1d(0,1)
        p1.x_range = bokeh.models.Range1d(-.1,7.1)
        if np.max(df_abundance_good_meso['relative_abundance'])*1.1>current_max:
            p1.y_range = bokeh.models.Range1d(0,np.max(df_abundance_good_meso['relative_abundance'])*1.1)
            
            current_max = np.max(df_abundance_good_meso['relative_abundance'])*1.1
      #  p1.xaxis.visible=False
        p2.legend.visible=False
        custom_ticks = [0, current_max]
        p1.xaxis.ticker = FixedTicker(ticks=custom_ticks)
        #p1.xaxis.ticker.desired_num_ticks = 5
        plots_abundance.append(p1)
        plots_strain.append(p2)
        type_mesos.append(good_meso)
    return plots_abundance,plots_strain,type_mesos
      #  bokeh.io.show(p)



def transform_df(df_abundance):
    df_abundance['Lineage'] = df_abundance['species_id'].transform(lambda x: df_metadata.loc[df_metadata['species_id'] == x,'Lineage'].values[0])
    df_abundance['species'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[-1])
    df_abundance['genus'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[-2])
    df_abundance['family'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[-3])
    df_abundance['phyla'] = df_abundance['Lineage'].transform(lambda x: x.split(';')[1])
    return df_abundance


def get_both_dfs(fname):
    e003_metadata = pd.read_csv('e003_coalescence_metadata_round4_good.csv').set_index('sample')
    df1 = pd.read_csv(fname).drop(columns='Unnamed: 0').set_index('sample')
    df2 = pd.read_csv(fname.split('parent1_info.csv')[0] + 'parent2_info.csv').drop(columns='Unnamed: 0').set_index('sample')
    df2['conf_int'] = df2['boot_high']-df2['boot_low']
    
    df1['conf_int'] = df1['boot_high']-df1['boot_low']
    df_both = pd.concat([df1.rename(columns = {'boot_med':'boot_med1', 'boot_low':'boot_low1', 'boot_high':'boot_high1',
       'actual_med':'actual_med1', 'conf_int':'conf_int1'}),
                         df2.rename(columns = {'boot_med':'boot_med2', 'boot_low':'boot_low2', 'boot_high':'boot_high2',
       'actual_med':'actual_med2', 'conf_int':'conf_int2'})],axis=1)
    
    df_both['species'] = fname.split('/')[-2]
    df_both['fname'] = fname
  #  print('-'.join(fname1.split('/')[-1].split('-')[:-1]))
   # parent_media = 
    
   # df_both['subjects_measured'] = '-'.join(fname1.split('/')[-1].split('-')[:-1])
    #df_both['in_measured'] = '-'.join(fname1.split('/')[-1].split('_')[:-2])
    df_both['total_shift'] = np.abs(1-(df_both['actual_med1'] + df_both['actual_med2']))
    df_both['total_shift12'] = np.abs(1-(df_both['boot_low1'] + df_both['boot_high2']))
    df_both['total_shift21'] = np.abs(1-(df_both['boot_low2'] + df_both['boot_high1']))
    df_both['total_shift_max'] = df_both['total_shift21']
    df_both.loc[df_both['total_shift_max']<df_both['total_shift12'],'total_shift_max'] = df_both.loc[df_both['total_shift_max']<df_both['total_shift12'],'total_shift12']

    return df_both


def plot_abundance_trajectories_with_in_abundance(freq_df, df_abundance, colorby='meso',alpha=1.,strain1=True,
                                                  shift=True,markers='meso',
                             p1=None,p2=None,emph=[]):
    if colorby=='media':
        cmap = {'fecal': bokeh.palettes.Set2[8][-2],
        'mBHI':bokeh.palettes.Set2[8][1],
        'mGAM': bokeh.palettes.Set2[8][0],}
    if colorby=='const':
        color = bokeh.palettes.Set2[8][2]
  #  marker_dic =['circle_cross','circle_dot','circle_x','circle_y','circle']
    marker_dic = ['circle','triangle','square', 'diamond', 'square_pin','triangle_pin']

    if not p2:
        
        p2 = bokeh.plotting.figure(width = 500, height = 300)
        p2.xaxis.axis_label_text_font_size=label_font_size
        p2.xaxis.major_label_text_font_size=tick_font_size
        p2.yaxis.axis_label_text_font_size=label_font_size
        p2.yaxis.major_label_text_font_size=tick_font_size
        p2.xaxis.minor_tick_line_color= None
    if not p1:
        p1 = bokeh.plotting.figure(width = 500, height = 150)
        p1.xaxis.axis_label_text_font_size=label_font_size
        p1.xaxis.major_label_text_font_size=tick_font_size
        p1.yaxis.axis_label_text_font_size=label_font_size
        p1.yaxis.major_label_text_font_size=tick_font_size
        p1.xaxis.minor_tick_line_color= None

    for i, type_meso in enumerate(freq_df['type_mesocosm'].unique()):
        mesos = freq_df.loc[freq_df['type_mesocosm'] ==type_meso, 'mesocosm'].unique()
        media = type_meso.split('-')[-1]
        if colorby=='media':
            color = cmap[media]
        if colorby=='meso':
            if media == 'mBHI':
                #cmap = [bokeh.palettes.Light[9][1], bokeh.palettes.Light[9][2], bokeh.palettes.Light[9][3], bokeh.palettes.Bright[7][2]]
               cmap = bokeh.palettes.Light[9][3:]
            else:
                cmap= [bokeh.palettes.Light[9][0], bokeh.palettes.Light[9][5], bokeh.palettes.Light[9][6], bokeh.palettes.Light[9][7]]
                cmap = bokeh.palettes.Light[9][3:]

        for j,mesocosm in enumerate(mesos):
            df_meso = freq_df.loc[freq_df['mesocosm'] == mesocosm,:]
            df_abundance_meso = df_abundance.loc[df_abundance['mesocosm']==mesocosm,:]
            in_sample = df_meso['inoculumn_sample'].unique()[0]
            rep = df_meso['replicate'].unique()[0]
            marker = 'circle'
            if colorby=='meso':
                color=cmap[int(rep)]
            if markers == 'meso':
                marker = marker_dic[int(rep)]

            if in_sample not in df_meso['sample'].values:
                df_meso= pd.concat([df_meso,freq_df.loc[freq_df['sample'] == in_sample,:]],axis=0)
                df_abundance_meso = pd.concat([df_abundance_meso, df_abundance.loc[df_abundance['sample']==in_sample]],axis=0)
            df_meso = df_meso.sort_values(by = 'passage') 
            df_abundance_meso = df_abundance_meso.sort_values(by='passage')
            if strain1:
                minor_strain = 'boot_med1'
                upper = 'boot_high1'
                lower = 'boot_low1'
                if shift:
                    df_meso['boot_med1_shift']=df_meso['boot_med1']/(df_meso['boot_med1']+df_meso['boot_med2'])
                    df_meso['boot_high1_shift']=df_meso['boot_high1']/(df_meso['boot_high1']+df_meso['boot_high2'])
                    df_meso['boot_low1_shift']=df_meso['boot_low1']/(df_meso['boot_low1']+df_meso['boot_low2'])
                    minor_strain = 'boot_med1_shift'
                    upper = 'boot_high1_shift'
                    lower = 'boot_low1_shift'
            else:
                minor_strain = 'boot_med2'
                upper = 'boot_high2'
                lower = 'boot_low2'
            df_source = bokeh.models.ColumnDataSource(data=dict(base=df_meso['passage'].values, 
                                                                upper=df_meso[upper].values, lower=df_meso[lower].values))
            error = bokeh.models.Whisker(base='base', upper='upper', lower='lower', source=df_source,line_color='black',
                                         line_width=1,  level="annotation",)
            error.upper_head.size=5
            error.lower_head.size=5
            p1.line(df_abundance_meso['passage'].values,
                   df_abundance_meso['relative_abundance'].values,color=color,alpha=alpha)
            p1.scatter(df_abundance_meso['passage'].values,
                   df_abundance_meso['relative_abundance'].values,color = color, legend_label = mesocosm, size = 8,
                   marker=marker,#line_color='black',
                      alpha=alpha)
            
            p2.line(df_meso['passage'].values,
                   df_meso[minor_strain].values,color=color,alpha=alpha)
            p2.scatter(df_meso['passage'].values,
                   df_meso[minor_strain].values,color = color, legend_label = mesocosm, size = 8,marker=marker,
                      alpha=alpha)
            p2.add_layout(error,)
    subs_dic = {'AA':'A','AE':'B','AF':'C'}
    sub1, sub2 = type_meso.split('-')[0],type_meso.split('-')[1]
    sub1 = subs_dic[sub1] 
    sub2 = subs_dic[sub2] 
    if strain1:
        p2.yaxis.axis_label = f'Freq Strain {sub1} vs {sub2}'
    else:
        p2.yaxis.axis_label = f'Freq Strain {sub2} vs {sub1}'

    p1.xaxis.axis_label = 'Timepoint'
    # p2.y_range = bokeh.models.Range1d(-.005,1.005)
    # p2.x_range = bokeh.models.Range1d(-.5,7.5)
    p1.xgrid.grid_line_color = None
    p1.ygrid.grid_line_color = None

    p1.yaxis.axis_label = 'Sp Rel Abun'
    p2.xaxis.axis_label = 'Timepoint'
    p1.xaxis.ticker.desired_num_ticks = 6
    # p2.y_range = bokeh.models.Range1d(-.005,1.005)
    # p2.x_range = bokeh.models.Range1d(-.5,7.5)
    p2.xgrid.grid_line_color = None
    p2.ygrid.grid_line_color = None

    return p1,p2

def make_plain_linear_plot_abundance(df_abundance, species, inoculumn,emph=[],same_plots = False,colorby='media',markers='meso',shift=True, medias=['mBHI','mGAM']):
    e003_metadata = pd.read_csv('e003_coalescence_metadata_round4_good.csv').set_index('sample')
   # print(e003_metadata.columns.values)

    fname1=f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/{inoculumn}_parent1_info.csv'
    df_both = get_both_dfs(fname1)
    df_both['boot_med1_shift'] = df_both['boot_med1']/(df_both['boot_med1']+df_both['boot_med2'])
    df_both = df_both.loc[np.intersect1d(df_both.index.values,e003_metadata.index.values),:]
    df_both_meta = pd.concat([df_both,e003_metadata.loc[np.intersect1d(df_both.index.values,
                                                                               e003_metadata.index.values),:]],axis=1).reset_index()
    df_both_meta_good = df_both_meta.loc[df_both_meta['total_shift']<.1,:]
    if shift:
        df_both_meta_good['boot_med1_shift']=df_both_meta_good['boot_med1']/(df_both_meta_good['boot_med1']+df_both_meta_good['boot_med2'])
        df_both_meta_good['boot_high1_shift']=df_both_meta_good['boot_high1']/(df_both_meta_good['boot_high1']+df_both_meta_good['boot_high2'])
        df_both_meta_good['boot_low1_shift']=df_both_meta_good['boot_low1']/(df_both_meta_good['boot_low1']+df_both_meta_good['boot_low2'])
    df_abundance_good = df_abundance.loc[df_abundance['species_id']==species,:]
    plots_abundance = []
    plots_strain=[]
    type_mesos=[]
    if same_plots:
        p2 = bokeh.plotting.figure(width = 500, height = 300)
        p2.xaxis.axis_label_text_font_size=label_font_size
        p2.xaxis.major_label_text_font_size=tick_font_size
        p2.yaxis.axis_label_text_font_size=label_font_size
        p2.yaxis.major_label_text_font_size=tick_font_size
        p2.xaxis.minor_tick_line_color= None
        p1 = bokeh.plotting.figure(width = 500, height = 150)
        p1.xaxis.axis_label_text_font_size=label_font_size
        p1.xaxis.major_label_text_font_size=tick_font_size
        p1.yaxis.axis_label_text_font_size=label_font_size
        p1.yaxis.major_label_text_font_size=tick_font_size
        p1.xaxis.minor_tick_line_color= None


    current_max = 0
    for good_meso in df_both_meta_good['type_mesocosm'].unique():
        med = good_meso.split('-')[-1]
        if med not in medias:
            continue
        df_both_meta_good_meso= df_both_meta_good.loc[(df_both_meta_good['type_mesocosm'] == good_meso)+(df_both_meta_good['is_inoculumn']),:].copy()
        df_abundance_good_meso = df_abundance_good.loc[(df_abundance_good['type_mesocosm'] == good_meso)+(df_abundance_good['is_inoculumn']),:].copy()
        df_abundance_good_meso['type_mesocosm'] = good_meso
        df_both_meta_good_meso['type_mesocosm'] = good_meso
        if same_plots:
            p1,p2 = plot_abundance_trajectories_with_in_abundance(df_both_meta_good_meso,df_abundance_good_meso, strain1=True,alpha=1,
                                     emph=emph,p1=p1,p2=p2,colorby=colorby,markers=markers,shift=shift)
        else:
            p1,p2 = plot_abundance_trajectories_with_in_abundance(df_both_meta_good_meso, df_abundance_good_meso, strain1=True,alpha=1,
                                     emph=emph,colorby=colorby,markers=markers,shift=shift)
    
        species = fname1.split('/')[-2]
       # species_name = df_metadata.loc[species,'species_plot']
        subject1 = fname1.split('/')[-1].split('-')[0]
        in_name = fname1.split('/')[-1].split('_')[0]
       # p2.title.text = f'{species_name}, {subject1} in {in_name} collisions'
        p2.x_range = bokeh.models.Range1d(-.1,7.1)
        p2.y_range = bokeh.models.Range1d(0,1)
        p1.x_range = bokeh.models.Range1d(-.1,7.1)
        p1.y_range = bokeh.models.Range1d(0,np.max(df_abundance_good_meso['relative_abundance'])*1.1)
        if same_plots:
            if np.max(df_abundance_good_meso['relative_abundance'])*1.1>current_max:
                p1.y_range = bokeh.models.Range1d(-.001,np.max(df_abundance_good_meso['relative_abundance'])*1.1)
                current_max = np.max(df_abundance_good_meso['relative_abundance'])*1.1
      #  p1.xaxis.visible=False
        p2.legend.visible=False
        plots_abundance.append(p1)
        plots_strain.append(p2)
        type_mesos.append(good_meso)
    return plots_abundance,plots_strain,type_mesos

def plot_trajectories_with_in(freq_df, colorby='const',alpha=1.,strain1=True,
                                                  shift=True,markers='markers',
                         p2=None,emph=[]):
    if colorby=='media':
        cmap = {'fecal': bokeh.palettes.Set2[8][-2],
        'mBHI':bokeh.palettes.Set2[8][1],
        'mGAM': bokeh.palettes.Set2[8][0],}
  #  marker_dic =['circle_cross','circle_dot','circle_x','circle_y','circle']
    marker_dic = ['circle','triangle','square', 'diamond', 'hexagon']
    if not p2:
        p2 = bokeh.plotting.figure(width = 500, height = 300)
        p2.xaxis.axis_label_text_font_size=label_font_size
        p2.xaxis.major_label_text_font_size=tick_font_size
        p2.yaxis.axis_label_text_font_size=label_font_size
        p2.yaxis.major_label_text_font_size=tick_font_size
        p2.xaxis.minor_tick_line_color= None

    for i, type_meso in enumerate(freq_df['type_mesocosm'].unique()):
        mesos = freq_df.loc[freq_df['type_mesocosm'] ==type_meso, 'mesocosm'].unique()
        media = type_meso.split('-')[-1]
        if colorby=='media':
            color = cmap[media]
        if colorby=='const':
            color = bokeh.palettes.Set2[8][2]
        if colorby=='meso':
            if media == 'mBHI':
                #cmap = [bokeh.palettes.Light[9][1], bokeh.palettes.Light[9][2], bokeh.palettes.Light[9][3], bokeh.palettes.Bright[7][2]]
               cmap = bokeh.palettes.Light[9][3:]
            else:
                cmap= [bokeh.palettes.Light[9][0], bokeh.palettes.Light[9][5], bokeh.palettes.Light[9][6], bokeh.palettes.Light[9][7]]
                cmap = bokeh.palettes.Light[9][3:]


        for j,mesocosm in enumerate(mesos):
            df_meso = freq_df.loc[freq_df['mesocosm'] == mesocosm,:]
            in_sample = df_meso['inoculumn_sample'].unique()[0]
            rep = df_meso['replicate'].unique()[0]
            marker = 'circle'
            if colorby=='meso':
                color=cmap[int(rep)]
            if markers == 'meso':
                marker = marker_dic[int(rep)]
            if in_sample not in df_meso['sample'].values:
                df_meso= pd.concat([df_meso,freq_df.loc[freq_df['sample'] == in_sample,:]],axis=0)
            df_meso = df_meso.sort_values(by = 'passage') 
            if strain1:
                minor_strain = 'boot_med1'
                upper = 'boot_high1'
                lower = 'boot_low1'
                if shift:
                    minor_strain = 'boot_med1_shift'
                    upper = 'boot_high1_shift'
                    lower = 'boot_low1_shift'
            else:
                minor_strain = 'boot_med2'
                upper = 'boot_high2'
                lower = 'boot_low2'
            df_source = bokeh.models.ColumnDataSource(data=dict(base=df_meso['passage'].values, 
                                                                upper=df_meso[upper].values, lower=df_meso[lower].values))
            error = bokeh.models.Whisker(base='base', upper='upper', lower='lower', source=df_source,line_color='black',
                                         line_width=1,  level="annotation",)
            error.upper_head.size=5
            error.lower_head.size=5
            p2.line(df_meso['passage'].values,
                   df_meso[minor_strain].values,color=color,alpha=alpha)
            p2.scatter(df_meso['passage'].values,
                   df_meso[minor_strain].values,color = color, legend_label = mesocosm, size = 8,marker=marker,
                      alpha=alpha)
            p2.add_layout(error,)
    subs_dic = {'AA':'A','AE':'B','AF':'C'}
    sub1, sub2 = type_meso.split('-')[0],type_meso.split('-')[1]
    sub1 = subs_dic[sub1] 
    sub2 = subs_dic[sub2] 
    if strain1:
        p2.yaxis.axis_label = f'Freq Strain {sub1} vs {sub2}'
    else:
        p2.yaxis.axis_label = f'Freq Strain {sub2} vs {sub1}'
    p2.xaxis.axis_label = 'Timepoint'
    p2.xgrid.grid_line_color = None
    p2.ygrid.grid_line_color = None

    p2.xaxis.axis_label_text_font_size=label_font_size
    p2.xaxis.major_label_text_font_size=tick_font_size
    p2.yaxis.axis_label_text_font_size=label_font_size
    p2.yaxis.major_label_text_font_size=tick_font_size
    p2.xaxis.minor_tick_line_color= None
    p2.yaxis.minor_tick_line_color= None
    
    return p2

def make_logit_plot(species, inoculumn,emph=[],same_plots = False,colorby='const',markers='meso',shift=True,thresh=1e-3,medias = ['mBHI','mGAM'],logit=True):
    e003_metadata = pd.read_csv('e003_coalescence_metadata_round4_good.csv').set_index('sample')
   # print(e003_metadata.columns.values)

    fname1=f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/{inoculumn}_parent1_info.csv'
    df_both = get_both_dfs(fname1)
    df_both = df_both.loc[np.intersect1d(df_both.index.values,e003_metadata.index.values),:]
    df_both_meta = pd.concat([df_both,e003_metadata.loc[np.intersect1d(df_both.index.values,
                                                                               e003_metadata.index.values),:]],axis=1).reset_index()

    df_both_meta_good = df_both_meta.loc[df_both_meta['total_shift']<.1,:]
    for col in ['boot_med1','boot_low1','boot_high1','boot_med2','boot_low2','boot_high2']:
        df_both_meta_good.loc[df_both_meta_good[col]<thresh,col]=thresh
        df_both_meta_good.loc[df_both_meta_good[col]>1-thresh,col]=1-thresh
    if shift:
        df_both_meta_good['boot_med1_shift']=df_both_meta_good['boot_med1']/(df_both_meta_good['boot_med1']+df_both_meta_good['boot_med2'])
        df_both_meta_good['boot_high1_shift']=df_both_meta_good['boot_high1']/(df_both_meta_good['boot_high1']+df_both_meta_good['boot_high2'])
        df_both_meta_good['boot_low1_shift']=df_both_meta_good['boot_low1']/(df_both_meta_good['boot_low1']+df_both_meta_good['boot_low2'])
        for col in ['boot_med1_shift','boot_low1_shift','boot_high1_shift']:
            df_both_meta_good.loc[df_both_meta_good[col]<thresh,col]=thresh
            df_both_meta_good.loc[df_both_meta_good[col]>1-thresh,col]=1-thresh

    sel_stuff, df_both_meta_good = get_selection_stuff_both(species,inoculumn)
    df_both_meta_good = df_both_meta.loc[df_both_meta['total_shift']<.1,:]
    df_both_meta_good_logit = df_both_meta_good.copy()
    if logit:
        df_both_meta_good_logit[['boot_med1','boot_low1','boot_high1',
        'boot_med2','boot_low2','boot_high2',]] = np.log(df_both_meta_good_logit[['boot_med1','boot_low1','boot_high1',
        'boot_med2','boot_low2','boot_high2',]]/(1-df_both_meta_good_logit[['boot_med1','boot_low1','boot_high1',
        'boot_med2','boot_low2','boot_high2',]]) )
       # if shift:
        #    df_both_meta_good_logit[['boot_med1_shift','boot_low1_shift','boot_high1_shift']] = \
         #   np.log(df_both_meta_good_logit[['boot_med1_shift','boot_low1_shift','boot_high1_shift']]/(1-\
          #   df_both_meta_good_logit[['boot_med1_shift','boot_low1_shift','boot_high1_shift']]))
    plots_strain=[]
    type_mesos=[]
    for good_meso in df_both_meta_good['type_mesocosm'].unique():
        med = good_meso.split('-')[-1]
        if med not in medias:
            continue
        df_both_meta_good_meso_logit= df_both_meta_good_logit.loc[(df_both_meta_good_logit['type_mesocosm'] == good_meso)+(df_both_meta_good_logit['is_inoculumn']),:].copy()
        df_both_meta_good_meso_logit['type_mesocosm'] = good_meso
        p= hv.Points(np.random.rand(50,2)).opts(width = 500, height = 300,color=None)
        if logit:
            exps_adjust = np.array([.001,.01,.1,.25, .5, .75, .9,.99,.999])
            lins_adjust = np.log(exps_adjust/(1-exps_adjust))
            p.opts(yticks=[(lins_adjust[i], exps_adjust[i]) for i in range(len(lins_adjust))]
                    )
        p=hv.render(p)

        p = plot_trajectories_with_in(df_both_meta_good_meso_logit, strain1=True,shift=False, colorby=colorby,p2=p)
        subject1 = fname1.split('/')[-1].split('-')[0]
        in_name = fname1.split('/')[-1].split('_')[0]
        p.x_range = bokeh.models.Range1d(-.1,7.1)
        p.y_range = bokeh.models.Range1d(0,1)
        if logit:
            exps_adjust = np.array([.001,.01,.1,.25,.5,.75,.9,.99,.999])
            p.yaxis.ticker = np.log(exps_adjust/(1-exps_adjust))
            p.y_range = bokeh.models.Range1d(np.log(thresh/(1-thresh)),np.log((1-thresh)/thresh))
        p.legend.visible=False
        plots_strain.append(p)
        type_mesos.append(good_meso)
        p.xaxis.axis_label_text_font_size=label_font_size
        p.xaxis.major_label_text_font_size=tick_font_size
        p.yaxis.axis_label_text_font_size=label_font_size
        p.yaxis.major_label_text_font_size=tick_font_size
        p.xaxis.minor_tick_line_color= None

    return plots_strain,type_mesos




def get_plotting_stuff(species,inoculumn, polarize=False, sel_period=(1,3),adjust_low_freq=True):
    fname1 = f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_sel_bootstrapv3/{species}/{inoculumn}_parent1_info.csv'
    fname2 = f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_sel_bootstrapv3/{species}/{inoculumn}_parent2_info.csv'
    sel_stuff1 = pd.read_csv(fname1).rename(columns={'mesocosms':'mesocosm'})
  #  sel_stuff2 = pd.read_csv(fname2).rename(columns={'mesocosms':'mesocosm'})
    
    fname1=f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/{inoculumn}_parent1_info.csv'
    fname2=f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/{inoculumn}_parent1_info.csv'
    freq1 = pd.read_csv(fname1)
    freq2=pd.read_csv(fname2)
    df_both = get_both_dfs(fname1)
    e003_metadata = pd.read_csv('e003_coalescence_metadata_round4_good.csv').set_index('sample')
    df_both = df_both.loc[np.intersect1d(df_both.index.values,e003_metadata.index.values),:]
    df_both=df_both.loc[df_both['total_shift']<.1,:]

    df_both_meta = pd.concat([df_both,e003_metadata.loc[np.intersect1d(df_both.index.values,
                                                                                   e003_metadata.index.values),:]],
                                         axis=1).reset_index()
    sel_stuff1=sel_stuff1.loc[sel_stuff1['sample1'].isin(df_both_meta['sample'].unique())*sel_stuff1['sample2'].isin(df_both_meta['sample'].unique()),:]
    sel_stuff1=sel_stuff1.loc[(sel_stuff1['passage1']==sel_period[0])*(sel_stuff1['passage2']==sel_period[1]),:]
   # sel_stuff2=sel_stuff2.loc[sel_stuff2['sample1'].isin(df_both_meta['sample'].unique())*sel_stuff2['sample2'].isin(df_both_meta['sample'].unique()),:]
    #sel_stuff2=sel_stuff2.loc[(sel_stuff2['passage1']==sel_period[0])*(sel_stuff2['passage2']==sel_period[1]),:]
    
    df_both_meta_pfinal = df_both_meta.loc[df_both_meta['passage']== 7,:].set_index('mesocosm')

    df_both_meta_p_pol = df_both_meta.loc[df_both_meta['passage']== 0,:].set_index('mesocosm')
    df_both_meta_p_polgr=df_both_meta_p_pol.groupby(['type_mesocosm']).median(numeric_only=True).reset_index()
    to_repol = df_both_meta_p_polgr['actual_med1']>.5
    if len(to_repol) <1:
        return pd.DataFrame()
#    print(to_repol)
    to_repol=to_repol[0]
  #  print(to_repol)
    
    sel_stuff1=sel_stuff1.set_index('mesocosm')
 #   sel_stuff2=sel_stuff2.set_index('mesocosm')

    good_mesos = np.intersect1d(sel_stuff1.index.values,df_both_meta_pfinal.index.values)
    sel_stuff1=sel_stuff1.loc[good_mesos,:]
    df_both_meta_pfinal=df_both_meta_pfinal.loc[good_mesos,:]
    df_both_meta_pfinal['sel_coeff']= np.nan
    df_both_meta_pfinal.loc[good_mesos,'self_coeff'] = sel_stuff1.loc[good_mesos,'sel_med']

    
    return df_both_meta_pfinal

    

def adjust_df(meso_df, good_families ):
    for i,fam in enumerate(good_families):
        meso_df.loc[(meso_df['family']==fam)*(meso_df['relative_abundance']<1e-2),'species_id'] = i
    meso_df = meso_df.groupby(['family','sample','passage','species_id',]).sum().reset_index()
                               
    for sample in meso_df['sample'].unique():
        full= meso_df.loc[meso_df['sample']==sample,'relative_abundance'].sum()
        if full>0:
            #print(sample,full)
            meso_df.loc[meso_df['sample']==sample,'relative_abundance']=  meso_df.loc[meso_df['sample']==sample,'relative_abundance']/full
    
    return meso_df

def make_bar_plot_assembly(meso_to_look_at, df_abundance,add=.2,bar_width = 1,plot_width=600,include_glycerol=False, species_emph = None):
    df_meso=df_abundance.loc[df_abundance['mesocosm']==meso_to_look_at,:]
    sub = meso_to_look_at.split('-')[1]
    in_sample = df_abundance.loc[df_abundance['mesocosm']==f'fecal-{sub}-fecal',:]
    df_meso = pd.concat([df_meso,in_sample])
    df_meso['passage_plot']=df_meso['passage'].astype(str)
    df_meso['family'] = df_meso['family'].transform(lambda x: x.split('f__')[-1])
    sums = df_meso.groupby(['sample']).sum(numeric_only=True)
    good_samples = sums.loc[sums['relative_abundance']>0.,:].index.values
    df_meso=df_meso.loc[df_meso['sample'].isin(good_samples),:]
    cmap_family= {'Bacteroidaceae': '#8dd3c7',
    'Enterobacteriaceae': '#ffffb3',
    'Porphyromonadaceae': '#bebada',
    'Peptoniphilaceae': '#fb8072',
    'Oscillospiraceae': '#80b1d3',
    'Enterococcaceae': '#fdb462',
    'Tannerellaceae': '#b3de69',
    'Lachnospiraceae': '#fccde5',
    'Veillonellaceae': '#bc80bd',
    'Peptostreptococcaceae': '#ccebc5',
    'Acidaminococcaceae': '#ffed6f',
    'other': '#d9d9d9'}
    good_families = list(cmap_family.keys())
    df_meso.loc[~df_meso['family'].isin(good_families),'family']='other'
    df_meso = adjust_df(df_meso, good_families)
    df_meso['family_plot']=df_meso['family'].copy()
    #df_meso.loc[~df_meso['family'].isin(good_families),'family_plot']='other'
    df_meso=df_meso.sort_values(by='family_plot')
    passages = [0,1,2,3,4,5]
    if include_glycerol:
        passages = [0,1,2,3,4,5,6]
    df_meso_small =df_meso.loc[df_meso['passage'].isin(passages),:].sort_values(by='family_plot')
    bars2 = hv.Bars(df_meso_small, kdims=[hv.Dimension('passage', values=passages), 'species_id',],
                vdims = ['relative_abundance','family_plot',])
    alpha = 1.
    alpha_bar = .5
    if species_emph:
        alpha=.1
        alpha_bar = .1
    

    bars2=bars2.opts(width=plot_width, height=365).opts(stacked=True,#alpha='relative_abundance',
                                        color='family_plot',
                                        cmap=cmap_family,
                                        alpha=alpha,
                                                ylim = (0,1),
                                                bar_width = bar_width,
                                        xlabel='Timepoint',
                                        ylabel='Relative Abundance',
                                            show_legend=True,legend_position='right',
                                            )#.sort(by='passage_plot')   
    sp_good = df_meso_small['species_id'].unique()
    df_meso_small_adjust = df_meso_small.copy()
    new_df = []
    for passage in df_meso_small['passage'].unique():
        df_meso_smallg = df_meso_small.loc[df_meso_small['passage']==passage,:].copy()
        df_meso_smallg['passage']=passage+add
        new_df.append(df_meso_smallg)
        df_meso_smallb = df_meso_small.loc[df_meso_small['passage']==passage,:].copy()
        df_meso_smallb['passage']=passage-add
        new_df.append(df_meso_smallb)
    new_df = pd.concat(new_df)
    df_meso_small_big = pd.concat([new_df, df_meso_small])
    to_overlay = []
    for sp in sp_good:
        if str(sp) == str(species_emph):
            to_overlay.append(hv.Area(df_meso_small_big.loc[df_meso_small_big['species_id']==sp,:].sort_values(by='passage'), 
            kdims=[hv.Dimension('passage',values=[0,1,2,3,4,5],)],
                                vdims=['relative_abundance','family_plot']).opts(alpha=1.,line_color='black',
                                                                                color=cmap_family[df_meso_small.loc[df_meso_small['species_id']==sp,'family_plot'].values[0]]))
        else:
            to_overlay.append(hv.Area(df_meso_small_big.loc[df_meso_small_big['species_id']==sp,:].sort_values(by='passage'), 
            kdims=[hv.Dimension('passage',values=[0,1,2,3,4,5],)],
                                vdims=['relative_abundance','family_plot']).opts(alpha=alpha_bar,line_color='grey',
                                                                                color=cmap_family[df_meso_small.loc[df_meso_small['species_id']==sp,'family_plot'].values[0]]))

    overlay = hv.Overlay(to_overlay)
    p=hv.Area.stack(overlay)

    if species_emph:
        p = hv.render((bars2*p).opts(xlim=(-0.25,5.5),ylim=(0,1)))
    else:
        p = hv.render((p*bars2).opts(xlim=(-0.25,5.5),ylim=(0,1)))
   # p.legend.visible = False
    p.xaxis.axis_label_text_font_size=label_font_size
    p.xaxis.major_label_text_font_size=tick_font_size
    p.yaxis.axis_label_text_font_size=label_font_size
    p.yaxis.major_label_text_font_size=tick_font_size
    p.xaxis.minor_tick_line_color= None
    p.yaxis.minor_tick_line_color= None
    return p


def make_bar_plot_coalescence(df_meso,add=.2,bar_width = 1, passages = [0,1,2,3,4,5,6,7], species_emph = None):
    df_meso['passage_plot']=df_meso['passage'].astype(str)
    df_meso['family'] = df_meso['family'].transform(lambda x: x.split('f__')[-1])
    sums = df_meso.groupby(['sample']).sum(numeric_only=True)
    good_samples = sums.loc[sums['relative_abundance']>0.,:].index.values
    df_meso=df_meso.loc[df_meso['sample'].isin(good_samples),:]
    cmap_family= {'Bacteroidaceae': '#8dd3c7',
    'Enterobacteriaceae': '#ffffb3',
    'Porphyromonadaceae': '#bebada',
    'Peptoniphilaceae': '#fb8072',
    'Oscillospiraceae': '#80b1d3',
    'Enterococcaceae': '#fdb462',
    'Tannerellaceae': '#b3de69',
    'Lachnospiraceae': '#fccde5',
    'Veillonellaceae': '#bc80bd',
    'Peptostreptococcaceae': '#ccebc5',
    'Acidaminococcaceae': '#ffed6f',
    'other': '#d9d9d9'}
    good_families = list(cmap_family.keys())
    df_meso.loc[~df_meso['family'].isin(good_families),'family']='other'
    df_meso = adjust_df(df_meso, good_families)
    df_meso['family_plot']=df_meso['family'].copy()
    #df_meso.loc[~df_meso['family'].isin(good_families),'family_plot']='other'
    df_meso=df_meso.sort_values(by='family_plot')
    df_meso_small =df_meso.loc[df_meso['passage'].isin(passages),:].sort_values(by='family_plot')
    bars2 = hv.Bars(df_meso_small, kdims=[hv.Dimension('passage', values=passages), 'species_id',],
                vdims = ['relative_abundance','family_plot',])
    alpha = 1
    alpha_bar=.5
    if species_emph:
        alpha = .1
        alpha_bar=.1

    bars2=bars2.opts(width=500, height=250,).opts(stacked=True,#alpha='relative_abundance',
                                        color='family_plot',
                                        cmap=cmap_family,
                                        alpha=alpha,
                                                ylim = (0,1),
                                                bar_width = bar_width,
                                        xlabel='Timepoint',
                                        ylabel='Relative Abundance',
                                            show_legend=True,legend_position='right',
                                            )#.sort(by='passage_plot')   
                                        
    sp_good = df_meso_small['species_id'].unique()
    df_meso_small_adjust = df_meso_small.copy()
    new_df = []
    for passage in df_meso_small['passage'].unique():
        df_meso_smallg = df_meso_small.loc[df_meso_small['passage']==passage,:].copy()
        df_meso_smallg['passage']=passage+add
        new_df.append(df_meso_smallg)
        df_meso_smallb = df_meso_small.loc[df_meso_small['passage']==passage,:].copy()
        df_meso_smallb['passage']=passage-add
        new_df.append(df_meso_smallb)
    new_df = pd.concat(new_df)
    df_meso_small_big = pd.concat([new_df, df_meso_small])
    to_overlay = []
    for sp in sp_good:
        if str(sp) == str(species_emph):
            to_overlay.append(hv.Area(df_meso_small_big.loc[df_meso_small_big['species_id']==sp,:].sort_values(by='passage'), 
            kdims=[hv.Dimension('passage',values=[0,1,2,3,4,5],)],
                                vdims=['relative_abundance','family_plot']).opts(alpha=1.,line_color='black',
                                                                                color=cmap_family[df_meso_small.loc[df_meso_small['species_id']==sp,'family_plot'].values[0]]))
        else:
            to_overlay.append(hv.Area(df_meso_small_big.loc[df_meso_small_big['species_id']==sp,:].sort_values(by='passage'), 
            kdims=[hv.Dimension('passage',values=[0,1,2,3,4,5],)],
                                vdims=['relative_abundance','family_plot']).opts(alpha=alpha_bar,line_color='grey',
                                                                                color=cmap_family[df_meso_small.loc[df_meso_small['species_id']==sp,'family_plot'].values[0]]))

    overlay = hv.Overlay(to_overlay)
    p=hv.Area.stack(overlay)
    if species_emph:
        p = hv.render((bars2*p).opts(xlim=(-0.25,7.5),ylim=(0,1)))
    else:
        p = hv.render((p*bars2).opts(xlim=(-0.25,7.5),ylim=(0,1)))
    p.legend.visible = False
    p.xaxis.axis_label_text_font_size=label_font_size
    p.xaxis.major_label_text_font_size=tick_font_size
    p.yaxis.axis_label_text_font_size=label_font_size
    p.yaxis.major_label_text_font_size=tick_font_size
    p.xaxis.minor_tick_line_color= None
    p.yaxis.minor_tick_line_color= None

    return p


goodsp = [102544,101433,100150]
bad = [101288,100057, 103188, 102281, 102645]
def get_selection_stuff_both(species,inoculumn,thresh=1e-3):
    fname1 = f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_sel_bootstrapv3/{species}/{inoculumn}_parentboth_info.csv'
    sel_stuff = pd.read_csv(fname1).rename(columns={'mesocosms':'mesocosm'})
    # workflow/report/track_snpsv2_sel_bootstrapv3/100910/AA-AE-mBHI_parentboth_info.csv
    fname1=f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/{inoculumn}_parent1_info.csv'
    e003_metadata = pd.read_csv('e003_coalescence_metadata_round4_good.csv').set_index('sample')
    freq1 = pd.read_csv(fname1)
    df_both = get_both_dfs(fname1)
    df_both = df_both.loc[np.intersect1d(df_both.index.values,e003_metadata.index.values),:]
    df_both=df_both.loc[df_both['total_shift']<.1,:]
    

    df_both_meta_good = pd.concat([df_both,e003_metadata.loc[np.intersect1d(df_both.index.values,
                                                                                   e003_metadata.index.values),:]], axis=1).reset_index()


    print('bef',np.max(df_both_meta_good['actual_med1']- df_both_meta_good['boot_med1']))
    df_both_meta_good['actual_med1']=df_both_meta_good['actual_med1']/(df_both_meta_good['actual_med1']+df_both_meta_good['actual_med2'])
    df_both_meta_good['boot_med1']=df_both_meta_good['boot_med1']/(df_both_meta_good['boot_med1']+df_both_meta_good['boot_med2'])
    df_both_meta_good['boot_high1']=df_both_meta_good['boot_high1']/(df_both_meta_good['boot_high1']+df_both_meta_good['boot_high2'])
    df_both_meta_good['boot_low1']=df_both_meta_good['boot_low1']/(df_both_meta_good['boot_low1']+df_both_meta_good['boot_low2'])
    for col in ['boot_med1','boot_low1','boot_high1','boot_med2','boot_low2','boot_high2']:
        df_both_meta_good.loc[df_both_meta_good[col]<thresh,col]=thresh
        df_both_meta_good.loc[df_both_meta_good[col]>1-thresh,col]=1-thresh

    print('aft',np.max(df_both_meta_good['actual_med1']- df_both_meta_good['boot_med1']))
    sel_stuff=sel_stuff.loc[sel_stuff['sample1'].isin(df_both_meta_good['sample'].unique())*sel_stuff['sample2'].isin(df_both_meta_good['sample'].unique()),:]
    sel_stuff['type_mesocosm'] = sel_stuff['mesocosm'].transform(lambda x: '-'.join(x.split('-')[1:]))
    sel_stuff['parent_subjects'] = sel_stuff['mesocosm'].transform(lambda x: '-'.join(x.split('-')[1:3]))
    sel_stuff['parent_media'] = sel_stuff['mesocosm'].transform(lambda x: x.split('-')[3])
    sel_stuff['media'] = sel_stuff['mesocosm'].transform(lambda x: x.split('-')[-1])

    return sel_stuff, df_both_meta_good


def get_pred_sel(ts,sel_coeff, start_value,):
    pred=np.exp(ts*sel_coeff)
    pred = pred*start_value/(1 - start_value + start_value*pred)
    return pred 


def make_logit_plot_with_preds(species, inoculumn,t0,t1):
    fname1=f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_ALL_bootstrapv3/{species}/{inoculumn}_parent1_info.csv'
    sel_stuff = f'/Users/sophiewalton/git/coalescence-pilot-mgx/workflow/report/track_snpsv2_sel_bootstrapv3/{species}/{inoculumn}_parent1_info.csv'
    df_both = get_both_dfs(fname1)
    sel_stuff =pd.read_csv(sel_stuff )
    df_both['boot_med1_shift'] = df_both['boot_med1']/(df_both['boot_med1']+df_both['boot_med2'])
    df_both = df_both.loc[np.intersect1d(df_both.index.values,e003_metadata.index.values),:]
    df_both_meta = pd.concat([df_both,e003_metadata.loc[np.intersect1d(df_both.index.values,
                                                                               e003_metadata.index.values),:]],
                                     axis=1).reset_index()
          #  print(len(df_both_meta))
    df_both_meta_good = df_both_meta.loc[df_both_meta['total_shift']<.1,:]
    sel_stuff=sel_stuff.loc[sel_stuff['sample1'].isin(df_both_meta_good['sample'].values)*sel_stuff['sample2'].isin(df_both_meta_good['sample'].values),:]
    plots=[]
    for good_meso in df_both_meta_good['type_mesocosm'].unique():
    #    print(good_meso)
        df_both_meta_good_meso= df_both_meta_good.loc[(df_both_meta_good['type_mesocosm'] == good_meso)+(df_both_meta_good['is_inoculumn']),:]
 
        
        sel_stuff['type_meso']=sel_stuff['mesocosms'].transform(lambda x: '-'.join(x.split('-')[1:]))
        
        sel_stuff_meso=sel_stuff.loc[sel_stuff['type_meso']==good_meso,:]
        sel_stuff_meso=sel_stuff_meso.set_index('sample2')
        df_both_meta_good_meso=df_both_meta_good_meso.set_index('sample')
     #   print(df_both_meta_good_meso.loc[sel_stuff_meso.index.values,'boot_med1'].reset_index())
        sel_stuff_meso['start_value']=df_both_meta_good_meso.loc[sel_stuff_meso.index.values,'boot_med1'] #.reset_index()
        sel_stuff_meso=sel_stuff_meso.reset_index()
        sel_stuff_meso=sel_stuff_meso.set_index('sample1')
        sel_stuff_meso['start_value1']=df_both_meta_good_meso.loc[sel_stuff_meso.index.values,'boot_med1']
        sel_stuff_meso=sel_stuff_meso.reset_index()
        df_both_meta_good_meso=df_both_meta_good_meso.reset_index()
        #print(sel_stuff_meso)
        if len(sel_stuff_meso)==0:
            continue
     #   print(sel_stuff)
        sel_stuff_meso_meds = sel_stuff_meso.groupby(['passage1','passage2','type_meso']).median(numeric_only=True).reset_index()
      
        sel_coeff = sel_stuff_meso_meds.loc[(sel_stuff_meso_meds['passage1']==t0)*(sel_stuff_meso_meds['passage2']==t1),'sel_med']
        if len(sel_coeff)==0:
            continue
        start_value = df_both_meta_good_meso.loc[df_both_meta_good_meso['passage']==t1,'boot_med1'].median()
        
        
        sel_coeffs = sel_stuff_meso.loc[(sel_stuff_meso['passage1']==t0)*(sel_stuff_meso['passage2']==t1),'sel_med'].values
      #  print(sel_stuff_meso.loc[(sel_stuff_meso['passage1']==t0)*(sel_stuff_meso['passage2']==t1),'sample1'].values)
        b=sel_stuff_meso.loc[(sel_stuff_meso['passage1']==t0)*(sel_stuff_meso['passage2']==t1),'sample2'].values
      #  print(df_both_meta_good_meso.loc[df_both_meta_good_meso['sample'].isin(b),'boot_med1'])
       # samples=sel_stuff_meso.loc[(sel_stuff_meso['passage1']==t0)*(sel_stuff_meso['passage2']==t1),'sample'].values
      #  print(sel_stuff_meso.loc[(sel_stuff_meso['passage1']==t0)*(sel_stuff_meso['passage2']==t1),['start_value1','start_value','sel_med']])
        start_values=sel_stuff_meso.loc[(sel_stuff_meso['passage1']==t0)*(sel_stuff_meso['passage2']==t1),'start_value'].values

        p = make_logit_plot(df_both_meta_good_meso, sel_coeff.values[0], 
                            start_value, t1,)
        ts=np.arange(9)
        pred_values = get_pred_sel(ts-t1,sel_coeff.values[0], start_value,)
        pred_values_logit = np.log(pred_values/(1-pred_values))
        p.line(ts, pred_values_logit,width=1.5, color = 'grey')

        pred_values_t7=get_pred_sel(7-t1,sel_coeffs, start_values,)
     #   max_7 = np.max(pred_values_t7)
        argmax_7 =np.argmax(pred_values_t7) 
        argmin_7 =np.argmin(pred_values_t7) 

        ts=np.arange(9)
        pred_valuesmax = get_pred_sel(ts-t1,sel_coeffs[argmax_7], start_values[argmax_7],)
        pred_values_logitmax = np.log(pred_valuesmax/(1-pred_valuesmax))

        pred_valuesmin = get_pred_sel(ts-t1,sel_coeffs[argmin_7], start_values[argmin_7],)
        pred_values_logitmin = np.log(pred_valuesmin/(1-pred_valuesmin))

        print(pred_valuesmax >pred_valuesmin)
        print(sel_coeffs[argmax_7], start_values[argmax_7],sel_coeffs[argmin_7], start_values[argmin_7],)
        p.line(ts,pred_values_logitmin,color='grey',line_dash='dashed')
        p.line(ts,pred_values_logitmax,color='grey',line_dash='dashed')

        species = fname1.split('/')[-2]
        species_name = df_metadata.loc[species,'species_plot']
        subject1 = fname1.split('/')[-1].split('-')[0]
        in_name = fname1.split('/')[-1].split('_')[0]
        p.title.text = f'{species_name}, {subject1} in {in_name} collisions'
        p.x_range = bokeh.models.Range1d(-.1,7.1)
        p.y_range = bokeh.models.Range1d(np.log(5e-3/(1-5e-3)),np.log((1-5e-3)/(5e-3)))
        p.legend.visible = False
        p.xaxis.axis_label_text_font_size=label_font_size
        p.xaxis.major_label_text_font_size=tick_font_size
        p.yaxis.axis_label_text_font_size=label_font_size
        p.yaxis.major_label_text_font_size=tick_font_size
        p.xaxis.minor_tick_line_color= None
        p.yaxis.minor_tick_line_color= None
      #  bokeh.io.show(p)
        plots.append(p)
    return plots
        