"""
cd /Users/lu/Documents/chang/psoralen/covariation
python ~/Documents/scripts/maf/covariation_tree.py 

To use this script, change the following variables as necessary.
starttree: ["ensembl_amniota23.tree", "hg38_100way.tree"]
specieslist: shown below
"""


import ete2 as ete
xistaltspeciestree = ["Homo_sapiens", "Pan_troglodytes", "Pongo_pygmaeus_abelii", "Nomascus_leucogenys", "Macaca_mulatta", "Macaca_fascicularis", "Papio_anubis", "Saimiri_boliviensis", "Otolemur_garnettii", "Tupaia_chinensis", "Spermophilus_tridecemlineatus", "Rattus_norvegicus", "Octodon_degus", "Oryctolagus_cuniculus", "Ochotona_princeps", "Vicugna_pacos", "Orcinus_orca", "Pteropus_alecto", "Myotis_davidii", "Sorex_araneus", "Loxodonta_africana", "Trichechus_manatus_latirostris", "Orycteropus_afer_afer"]
xistspecieslist = ["Homo_sapiens", "Gorilla_gorilla_gorilla", "Pongo_pygmaeus_abelii", "Nomascus_leucogenys", "Macaca_mulatta", "Macaca_fascicularis", "Papio_anubis", "Chlorocebus_sabaeus", "Callithrix_jacchus", "Saimiri_boliviensis", "Otolemur_garnettii", "Tupaia_chinensis", "Spermophilus_tridecemlineatus", "Jaculus_jaculus", "Microtus_ochrogaster", "Cricetulus_griseus", "Mesocricetus_auratus", "Mus_musculus", "Rattus_norvegicus", "Heterocephalus_glaber", "Cavia_porcellus", "Chinchilla_lanigera", "Octodon_degus", "Oryctolagus_cuniculus", "Ochotona_princeps", "Sus_scrofa", "Vicugna_pacos", "Camelus_ferus", "Tursiops_truncatus", "Orcinus_orca", "Pantholops_hodgsonii", "Bos_taurus", "Ovis_aries", "Capra_hircus", "Ceratotherium_simum", "Felis_catus", "Canis_lupus_familiaris", "Mustela_putorius_furo", "Ailuropoda_melanoleuca", "Odobenus_rosmarus_divergens", "Leptonychotes_weddellii", "Pteropus_alecto", "Pteropus_vampyrus", "Eptesicus_fuscus", "Myotis_lucifugus", "Sorex_araneus", "Condylura_cristata", "Elephantulus_edwardii", "Trichechus_manatus_latirostris", "Chrysochloris_asiatica", "Orycteropus_afer_afer", "Dasypus_novemcinctus"]
cd9specieslist = ["Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla_gorilla", "Pongo_pygmaeus_abelii", "Nomascus_leucogenys", "Macaca_mulatta", "Macaca_fascicularis", "Papio_anubis", "Chlorocebus_sabaeus", "Callithrix_jacchus", "Saimiri_boliviensis", "Otolemur_garnettii", "Tupaia_chinensis", "Spermophilus_tridecemlineatus", "Jaculus_jaculus", "Microtus_ochrogaster", "Cricetulus_griseus", "Mus_musculus", "Rattus_norvegicus", "Heterocephalus_glaber", "Cavia_porcellus", "Chinchilla_lanigera", "Octodon_degus", "Oryctolagus_cuniculus", "Ochotona_princeps", "Sus_scrofa", "Vicugna_pacos", "Camelus_ferus", "Tursiops_truncatus", "Orcinus_orca", "Pantholops_hodgsonii", "Bos_taurus", "Ovis_aries", "Capra_hircus", "Equus_caballus", "Ceratotherium_simum", "Felis_catus", "Canis_lupus_familiaris", "Mustela_putorius_furo", "Ailuropoda_melanoleuca", "Odobenus_rosmarus_divergens", "Leptonychotes_weddellii", "Pteropus_alecto", "Pteropus_vampyrus", "Eptesicus_fuscus", "Myotis_davidii", "Myotis_lucifugus", "Erinaceus_europaeus", "Sorex_araneus", "Condylura_cristata", "Loxodonta_africana", "Elephantulus_edwardii", "Trichechus_manatus_latirostris", "Chrysochloris_asiatica", "Echinops_telfairi", "Dasypus_novemcinctus", "Monodelphis_domestica", "Sarcophilus_harrisii", "Macropus_eugenii", "Ornithorhynchus_anatinus", "Columba_livia", "Falco_cherrug", "Falco_peregrinus", "Ficedula_albicollis", "Geospiza_fortis", "Taeniopygia_guttata", "Pseudopodoces_humilis", "Melopsittacus_undulatus", "Amazona_vittata", "Ara_macao", "Anas_platyrhynchos", "Gallus_gallus", "Meleagris_gallopavo", "Alligator_mississippiensis", "Chelonia_mydas", "Chrysemys_picta_bellii", "Pelodiscus_sinensis", "Apalone_spinifera", "Xenopus_tropicalis", "Latimeria_chalumnae"]
krt5specieslist = ["Homo_sapiens","Pan_troglodytes","Gorilla_gorilla_gorilla","Pongo_pygmaeus_abelii","Nomascus_leucogenys","Macaca_mulatta","Macaca_fascicularis","Papio_anubis","Chlorocebus_sabaeus","Callithrix_jacchus","Saimiri_boliviensis","Otolemur_garnettii","Tupaia_chinensis","Spermophilus_tridecemlineatus","Jaculus_jaculus","Microtus_ochrogaster","Mesocricetus_auratus","Mus_musculus","Rattus_norvegicus","Heterocephalus_glaber","Cavia_porcellus","Chinchilla_lanigera","Octodon_degus","Oryctolagus_cuniculus","Sus_scrofa","Vicugna_pacos","Camelus_ferus","Orcinus_orca","Pantholops_hodgsonii","Bos_taurus","Ovis_aries","Capra_hircus","Equus_caballus","Ceratotherium_simum","Felis_catus","Canis_lupus_familiaris","Mustela_putorius_furo","Ailuropoda_melanoleuca","Odobenus_rosmarus_divergens","Leptonychotes_weddellii","Pteropus_alecto","Pteropus_vampyrus","Eptesicus_fuscus","Myotis_davidii","Myotis_lucifugus","Loxodonta_africana","Elephantulus_edwardii","Trichechus_manatus_latirostris","Chrysochloris_asiatica","Echinops_telfairi","Orycteropus_afer_afer","Dasypus_novemcinctus","Sarcophilus_harrisii"]
rpl8specieslist = ["Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla_gorilla", "Pongo_pygmaeus_abelii", "Nomascus_leucogenys", "Macaca_mulatta", "Macaca_fascicularis", "Papio_anubis", "Chlorocebus_sabaeus", "Callithrix_jacchus", "Otolemur_garnettii", "Tupaia_chinensis", "Spermophilus_tridecemlineatus", "Microtus_ochrogaster", "Cricetulus_griseus", "Heterocephalus_glaber", "Cavia_porcellus", "Chinchilla_lanigera", "Oryctolagus_cuniculus", "Sus_scrofa", "Vicugna_pacos", "Camelus_ferus", "Tursiops_truncatus", "Orcinus_orca", "Pantholops_hodgsonii", "Bos_taurus", "Ovis_aries", "Capra_hircus", "Equus_caballus", "Ceratotherium_simum", "Felis_catus", "Canis_lupus_familiaris", "Mustela_putorius_furo", "Ailuropoda_melanoleuca", "Odobenus_rosmarus_divergens", "Leptonychotes_weddellii", "Pteropus_alecto", "Pteropus_vampyrus", "Eptesicus_fuscus", "Myotis_davidii", "Myotis_lucifugus", "Loxodonta_africana", "Trichechus_manatus_latirostris", "Chrysochloris_asiatica", "Orycteropus_afer_afer", "Dasypus_novemcinctus", "Monodelphis_domestica", "Sarcophilus_harrisii", "Macropus_eugenii", "Ornithorhynchus_anatinus", "Columba_livia", "Falco_cherrug", "Falco_peregrinus", "Zonotrichia_albicollis", "Geospiza_fortis", "Taeniopygia_guttata", "Pseudopodoces_humilis", "Melopsittacus_undulatus", "Amazona_vittata", "Ara_macao", "Anas_platyrhynchos", "Gallus_gallus", "Meleagris_gallopavo", "Alligator_mississippiensis", "Chelonia_mydas", "Chrysemys_picta_bellii", "Pelodiscus_sinensis", "Apalone_spinifera", "Anolis_carolinensis", "Latimeria_chalumnae", "Tetraodon_nigroviridis", "Takifugu_rubripes", "Takifugu_flavidus", "Oreochromis_niloticus", "Neolamprologus_brichardi", "Maylandia_zebra", "Pundamilia_nyererei", "Oryzias_latipes", "Xiphophorus_maculatus", "Gasterosteus_aculeatus", "Gadus_morhua", "Astyanax_mexicanus", "Lepisosteus_oculatus", "Petromyzon_marinus"]
specieslist = xistaltspeciestree
#starttree = "ensembl_amniota23.tree"
starttree = "hg38_100way.tree"
outpdf = "hsXIST_alt1_tree.pdf"




tree = ete.Tree(starttree)
tree.prune(specieslist, preserve_branch_length = False)

for n in tree.traverse():
    style = ete.NodeStyle()
    #style['hz_line_width'] = 1
    #style['vt_line_width'] = 1
    style['size'] = 0
    n.set_style(style)
ts = ete.TreeStyle()
ts.mode = 'r'
ts.show_leaf_name = False
ts.show_scale = False
tree.render(outpdf, tree_style = ts, h=200)
tree.show(tree_style = ts)




"""
grep -f <(grep ">" hg38_multiz100_RPL8_group1_mlocarna.fa | \
cut -c2-) ~/Documents/chang/psoralen/covariation/hg38_100way.name | \
cut -f5 | tr '\n' ',' | sed 's/,/", "/g'

multiz100tree.prune(rpl8specieslist, preserve_branch_length = False)

for n in multiz100tree.traverse():
    style = ete.NodeStyle()
    #style['hz_line_width'] = 1
    #style['vt_line_width'] = 1
    style['size'] = 0
    n.set_style(style)

ts = ete.TreeStyle()
ts.mode = 'r'
ts.show_leaf_name = False
ts.show_scale = False

multiz100tree.render("RPL8_multiz100_tree.pdf", tree_style = ts, h=200) #, , w=50, units="mm")
multiz100tree.show(tree_style=ts)

"""










"""
example scripts from this website:
http://www.r-bloggers.com/phylogenies-in-r-and-python/

# load data
#traits = pd.read_csv('/Users/Nate/Documents/FIU/Research/Invasion_TraitPhylo/Data/plantTraits.csv')

#### TRAIT CLEANUP ####
# put an underscore in trait species
traits['species'] = traits['species'].map(lambda x: x.replace(' ', '_'))
# pull out the relevant traits and only keep complete cases
traits = traits[['species', 'introduced', 'woody', 'SLA', 'seedMass', 'toughness']]
traits = traits.dropna()
 
# next, prune down the traits data
traitsPrune = traits[traits['species'].isin(SERCphylo.get_leaf_names())]
 
# prune the phylogeny so only species with traits are kept
SERCphylo.prune(traitsPrune['species'], preserve_branch_length = True)


# guide for color
cols = [['black', 'red'][x] for x in traitsPrune['introduced']]
colorGuide = dict(zip(traitsPrune['species'], cols))
# weights (scaled to 1)
slaGuide = dict(zip(traitsPrune['species'], traitsPrune['SLA']/traitsPrune['SLA'].max()))
toughGuide = dict(zip(traitsPrune['species'], traitsPrune['toughness']/traitsPrune['toughness'].max()))
seedGuide = dict(zip(traitsPrune['species'], traitsPrune['seedMass']/traitsPrune['seedMass'].max()))

# set the base style of the phylogeny with thick lines
for n in SERCphylo.traverse():
style = ete.NodeStyle()
style['hz_line_width'] = 2
style['vt_line_width'] = 2
style['size'] = 0
n.set_style(style)

def mylayout(node):
# If node is a leaf, split the name and paste it back together to remove the underscore
if node.is_leaf():
temp = node.name.split('_')
sp = temp[0] + ' ' + temp[1]
temp2 = ete.faces.TextFace(sp, fgcolor = colorGuide[node.name], fsize = 18, fstyle = 'italic')


ts = ete.TreeStyle()
ts.mode = 'r'
ts.show_leaf_name = False
ts.layout_fn = mylayout
ts.branch_vertical_margin = 4
#ts.force_topology = True
ts.show_scale = False

def mylayout(node):
# If node is a leaf, split the name and paste it back together to remove the underscore
if node.is_leaf():
# species name
temp = node.name.split('_')
sp = temp[0] + ' ' + temp[1]
temp2 = ete.faces.TextFace(sp, fgcolor = colorGuide[node.name], fsize = 18, fstyle = 'italic')
ete.faces.add_face_to_node(temp2, node, column=0)
# make a circle for SLA, weighted by SLA values
sla = ete.CircleFace(radius = slaGuide[node.name]*15, color = colorGuide[node.name], style = 'circle')
sla.margin_left = 10
sla.hz_align = 1
ete.faces.add_face_to_node(sla, node, column = 0, position = 'aligned')
# same with toughness
toughness = ete.CircleFace(radius = toughGuide[node.name]*15, color = colorGuide[node.name], style = 'circle')
toughness.margin_left = 40
toughness.hz_align = 1
ete.faces.add_face_to_node(toughness, node, column = 1, position = 'aligned')

#/Users/zhipeng/Documents/chang/psoralen/covariation/hg38_100way_scientificnames.tree
"""
