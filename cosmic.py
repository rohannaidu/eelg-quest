
# coding: utf-8

# In[1]:

import numpy as np
import scipy
from scipy import stats
from astropy.io import fits
cosmos_list = fits.open('cosmos_3dhst.v4.1.cat.FITS', memmap=True)
cosmos_data = cosmos_list[1].data
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
pha_list = fits.open('3dhst.v4.1.4.full.v1.fits', memmap=True)
pha_data = pha_list[1].data
import numpy.ma as ma

def magab_err(flux, error):
    if (np.absolute(error)) >= flux:
        return 1.085
    else:
        return np.absolute((-1.085)* (error/flux))
    
def magab(flux, error):
    if error > flux:
        return (((-2.5)*(np.log10(error))) + 25 )
    else:
        return (((-2.5)*(np.log10(flux))) + 25 )

magfunc = np.vectorize(magab)
magerrorfunc = np.vectorize(magab_err)


# In[3]:

print scipy.stats.itemfreq(pha_data['id'][41200:75079] == cosmos_data['id'])
cosmos_phot_use_1= (cosmos_data['use_phot']==1)
print scipy.stats.itemfreq((cosmos_data['use_phot']==1))


# In[4]:

def cutter (a, a_err, b, b_err, c, c_err): 
    return (((magfunc(a, a_err)- magfunc(b, b_err) >(0.44 + (((magerrorfunc(a, a_err))**2) + ((magerrorfunc(b, b_err))**2))**.5)) & 
(magfunc(c, c_err)- magfunc(b, b_err) >(0.44 + (((magerrorfunc(c, c_err))**2) + ((magerrorfunc(b, b_err))**2))**.5))) & cosmos_phot_use_1 & ((a>=0)& (b>=0) & (c>=0)) & (((a/a_err)>3)&((b/b_err)>3)&((c/c_err)>3)))
#& (((a/a_err)>3)&((b/b_err)>3)&((c/c_err)>3))


# In[5]:

###jhk selection###
afils = [cosmos_data['f_j1'], cosmos_data['e_j1'], cosmos_data['f_j2'], cosmos_data['e_j2'], cosmos_data['f_j3'], cosmos_data['e_j3'] , cosmos_data['f_f125w'],cosmos_data['e_f125w'], cosmos_data['f_j'], cosmos_data['e_j'], cosmos_data['f_UVISTA_J'], cosmos_data['e_UVISTA_J'] ]
bfils= [cosmos_data['f_f140w'], cosmos_data['e_f140w'], cosmos_data['f_f160w'], cosmos_data['e_f160w'], cosmos_data['f_h1'], cosmos_data['e_h1'] , cosmos_data['f_h2'],cosmos_data['e_h2'], cosmos_data['f_h'], cosmos_data['e_h'], cosmos_data['f_UVISTA_H'], cosmos_data['e_UVISTA_H'] ]
cfils= [cosmos_data['f_k'], cosmos_data['e_k'], cosmos_data['f_ks'], cosmos_data['e_ks'], cosmos_data['f_uvista_ks'], cosmos_data['e_uvista_ks']]

for a in range(0,6):
    for b in range (0,6):
        for c in range (0,3):
            globals()['cosmos_jhk%d%d%d' %(a,b,c)] = (cutter(afils[2*a], afils[1+(2*a)], bfils[2*b], bfils[1+(2*b)], cfils[2*c], cfils[1+(2*c)]))

for a in range(0,6):
    for b in range (0,6):
        for c in range (0,3):
            print "cosmos_jhk%d%d%d"%(a,b,c), np.count_nonzero((globals()['cosmos_jhk%d%d%d' %(a,b,c)])==1)


# In[7]:

print (cosmos_data['id'][((cosmos_jhk302) | (cosmos_jhk311) | (cosmos_jhk312))])


# In[8]:

cosmos_oxy = pha_data['OIII_EQW'][41200:75079]
cosmos_hbeta = pha_data['HBETA_EQW'][41200:75079]
cosmos_eqw_sum = cosmos_oxy + cosmos_hbeta
cosmos_eqw_err_1 = (pha_data['OIII_EQW_ERR'][41200:75079])**2
cosmos_eqw_err_2 = (pha_data['HBETA_EQW_ERR'][41200:75079])**2
cosmos_eqw_err_net = (cosmos_eqw_err_1 + cosmos_eqw_err_2)**0.5

current_cut = (cosmos_jhk310|cosmos_jhk311|cosmos_jhk312)
#insert cut above
plt.figure(figsize=(16,8))
n_groups = len(cosmos_eqw_sum[(current_cut == 1)&((pha_data['HALPHA_EQW'][41200:75079] >0) | (cosmos_eqw_sum>0))])
print n_groups
print np.count_nonzero(((current_cut == 1)& ((pha_data['HALPHA_EQW'][41200:75079] >500) | (cosmos_eqw_sum>500)))==1)
print np.count_nonzero(((current_cut == 1)& ((pha_data['HALPHA_EQW'][41200:75079] >1000) | (cosmos_eqw_sum>1000)))==1)
means_men = cosmos_eqw_sum[(current_cut == 1)& ((pha_data['HALPHA_EQW'][41200:75079] >0) | (cosmos_eqw_sum>0))]
std_men = cosmos_eqw_err_net[(current_cut == 1)& ((pha_data['HALPHA_EQW'][41200:75079] >0) | (cosmos_eqw_sum>0))]
means_women = pha_data['HALPHA_EQW'][41200:75079][(current_cut == 1)& ((pha_data['HALPHA_EQW'][41200:75079] >0) | (cosmos_eqw_sum>0))]
std_women = pha_data['HALPHA_EQW_ERR'][41200:75079][(current_cut == 1) & ((pha_data['HALPHA_EQW'][41200:75079] >0) | (cosmos_eqw_sum>0))]
index = np.arange(n_groups)
bar_width = 0.35
opacity = 0.4
error_config = {'ecolor': '0.3'}
rects1 = plt.bar(index, means_men, bar_width,
                 alpha=opacity,
                 color='r',
                 yerr=std_men,
                 error_kw=error_config,
                 label='$OIII + H_{beta}$')
rects2 = plt.bar(index + bar_width, means_women, bar_width,
                 alpha=opacity,
                 color='b',
                 yerr=std_women,
                 error_kw=error_config,
                 label='$H_{alpha}$')
plt.xlabel('$ID$', fontsize = 25 )
plt.ylabel('$Equivalent Width$', fontsize = 25)
plt.title('$Cosmos JHK_{310|311|312} $', fontsize = 25)
plt.xticks(index + bar_width, pha_data['id'][41200:75079][(current_cut == 1)&((pha_data['HALPHA_EQW'][41200:75079] >0) | (cosmos_eqw_sum>0))])
plt.yticks(np.arange(0, max(means_men), 500 ))
plt.legend(loc='upper left')
plt.axis('tight')
plt.ylim(0,2000)


# In[9]:

scipy.stats.itemfreq(pha_data['field'][(((pha_data['OIII_EQW'] + pha_data['HBETA_EQW']) > 500) | (pha_data['HALPHA_EQW'] > 500))])


# In[10]:

scipy.stats.itemfreq(pha_data['field'][(((pha_data['z_spec']>1.75)&(pha_data['z_spec']<2.65)) | ((pha_data['z_peak']>1.75)&(pha_data['z_peak']<2.6))) & (((pha_data['OIII_EQW'] + pha_data['HBETA_EQW']) > 500) | (pha_data['HALPHA_EQW'] > 500))])


# In[4]:

if 1==1:
    for a in range(0,6):
        for b in range (0,6):
            for c in range (0,3):
                f1,f1e, f2, f2e, f3, f3e = afils[2*a], afils[1+(2*a)], bfils[2*b], bfils[1+(2*b)], cfils[2*c], cfils[1+(2*c)]
                plt.figure(figsize=(10,6))
                matplotlib.pyplot.scatter((magfunc(f2, f2e)- magfunc(f3, f3e)), (magfunc(f1, f1e)- magfunc(f2, f2e)) , 
                              c=['DarkOrange' if x==1&y==1 else 'r' if x==1 else 'LawnGreen' if y==1 else 'k' for x, y in zip(globals()['cosmos_jhk%d%d%d' %(a,b,c)],((((pha_data['z_spec'[41200:75079]]>1.75)&(pha_data['z_spec'][41200:75079]<2.65)) | ((pha_data['z_peak'][41200:75079]>1.75)&(pha_data['z_peak'][41200:75079]<2.6))) & (((pha_data['OIII_EQW'][41200:75079] + pha_data['HBETA_EQW'][41200:75079]) > 500) | (pha_data['HALPHA_EQW'][41200:75079] > 500))))], 
                              s= [50 if x==1&y==1 else 50 if x==1 else 50 if y==1 else 0.1 for x, y in zip(globals()['cosmos_jhk%d%d%d' %(a,b,c)],((((pha_data['z_spec'][41200:75079]>1.75)&(pha_data['z_spec'][41200:75079]<2.65)) | ((pha_data['z_peak'][41200:75079]>1.75)&(pha_data['z_peak'][41200:75079]<2.6))) & (((pha_data['OIII_EQW'][41200:75079] + pha_data['HBETA_EQW'][41200:75079]) > 500) | (pha_data['HALPHA_EQW'][41200:75079] > 500))))])
                plt.xlim(-5,5)
                plt.ylim(-1,5)
                plt.xlabel("$cosmos_jhk%d%d%d$"%(a,b,c), fontsize=16)
                plt.title('$Cosmos JHK$', fontsize=20)
                plt.savefig('cosmos_jhk%d%d%d.png'%(a,b,c),bbox_inches='tight')
                #save the plot to some folder rather than output
                # ssh turtle.astro.yale.edu
                # ssh terrapin.astro.yale.edu


# In[ ]:




# In[12]:

###ijh selection###
bfils = [cosmos_data['f_j1'], cosmos_data['e_j1'], cosmos_data['f_j2'], cosmos_data['e_j2'], cosmos_data['f_j3'], cosmos_data['e_j3'] , cosmos_data['f_f125w'],cosmos_data['e_f125w'], cosmos_data['f_j'], cosmos_data['e_j'], cosmos_data['f_UVISTA_J'], cosmos_data['e_UVISTA_J'] ]
cfils= [cosmos_data['f_f140w'], cosmos_data['e_f140w'], cosmos_data['f_f160w'], cosmos_data['e_f160w'], cosmos_data['f_h1'], cosmos_data['e_h1'] , cosmos_data['f_h2'],cosmos_data['e_h2'], cosmos_data['f_h'], cosmos_data['e_h'], cosmos_data['f_UVISTA_H'], cosmos_data['e_UVISTA_H'] ]
afils= [cosmos_data['f_I'], cosmos_data['e_I'], cosmos_data['f_Ip'], cosmos_data['e_Ip'], cosmos_data['f_f814W'], cosmos_data['e_f814W']]

for a in range(0,3):
    for b in range (0,6):
        for c in range (0,6):
            globals()['cosmos_ijh%d%d%d' %(a,b,c)] = (cutter(afils[2*a], afils[1+(2*a)], bfils[2*b], bfils[1+(2*b)], cfils[2*c], cfils[1+(2*c)]))

for a in range(0,3):
    for b in range (0,6):
        for c in range (0,6):
            print "cosmos_ijh%d%d%d"%(a,b,c), np.count_nonzero((globals()['cosmos_ijh%d%d%d' %(a,b,c)])==1)


# In[16]:

print (cosmos_data['id'][cosmos_ijh231 ])


# In[17]:

###ihk selection###
cfils = [cosmos_data['f_k'], cosmos_data['e_k'], cosmos_data['f_ks'], cosmos_data['e_ks'], cosmos_data['f_uvista_ks'], cosmos_data['e_uvista_ks']]
bfils= [cosmos_data['f_f140w'], cosmos_data['e_f140w'], cosmos_data['f_f160w'], cosmos_data['e_f160w'], cosmos_data['f_h1'], cosmos_data['e_h1'] , cosmos_data['f_h2'],cosmos_data['e_h2'], cosmos_data['f_h'], cosmos_data['e_h'], cosmos_data['f_UVISTA_H'], cosmos_data['e_UVISTA_H'] ]
afils= [cosmos_data['f_I'], cosmos_data['e_I'], cosmos_data['f_Ip'], cosmos_data['e_Ip'], cosmos_data['f_f814W'], cosmos_data['e_f814W']]

for a in range(0,3):
    for b in range (0,6):
        for c in range (0,3):
            globals()['cosmos_ihk%d%d%d' %(a,b,c)] = (cutter(afils[2*a], afils[1+(2*a)], bfils[2*b], bfils[1+(2*b)], cfils[2*c], cfils[1+(2*c)]))

for a in range(0,3):
    for b in range (0,6):
        for c in range (0,3):
            print "cosmos_ihk%d%d%d"%(a,b,c), np.count_nonzero((globals()['cosmos_ihk%d%d%d' %(a,b,c)])==1)


# In[18]:

print (cosmos_data['id'][cosmos_ihk210|cosmos_ihk211|cosmos_ihk212 ])


# In[19]:

###vij selection###
cfils = [cosmos_data['f_j1'], cosmos_data['e_j1'], cosmos_data['f_j2'], cosmos_data['e_j2'], cosmos_data['f_j3'], cosmos_data['e_j3'] , cosmos_data['f_f125w'],cosmos_data['e_f125w'], cosmos_data['f_j'], cosmos_data['e_j'], cosmos_data['f_UVISTA_J'], cosmos_data['e_UVISTA_J'] ]
bfils= [cosmos_data['f_I'], cosmos_data['e_I'], cosmos_data['f_Ip'], cosmos_data['e_Ip'], cosmos_data['f_f814W'], cosmos_data['e_f814W']]
afils= [cosmos_data['f_v'], cosmos_data['e_v'], cosmos_data['f_f606w'], cosmos_data['e_f606w']]

for a in range(0,2):
    for b in range (0,3):
        for c in range (0,6):
            globals()['cosmos_vij%d%d%d' %(a,b,c)] = (cutter(afils[2*a], afils[1+(2*a)], bfils[2*b], bfils[1+(2*b)], cfils[2*c], cfils[1+(2*c)]))

for a in range(0,2):
    for b in range (0,3):
        for c in range (0,5):
            print "cosmos_vij%d%d%d"%(a,b,c), np.count_nonzero((globals()['cosmos_vij%d%d%d' %(a,b,c)])==1)


# In[20]:

print (cosmos_data['id'][cosmos_vij123])


# In[ ]:



