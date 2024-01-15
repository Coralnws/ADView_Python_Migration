'''
Parse_Newick add to RT and TC class (alter some detail)
Read reference tree and tree collection separately (parsor will be a bit difference)
Result:
1. should have a reference tree with tc and all tc belongs to it will point to the reference tree
'''

from decimal import Decimal, getcontext
from Utils import *
getcontext().prec = 20

class Node:
    x = None # X-Coordinate
    y = None # Y-Coordinate

    def __init__(self, name, length=0):
        self.name = name
        self.length = length
        self.children = []
        self.support = None
        self.type = type  # Leaf or internal node

class Tree:
    root = None
    name = None
    node_list = []

    def __init__(self,name=None):
        self.name = name

    def print_tree(self, node=None, level=0):
        if node is None:
            node = self.root
        if node.name:
            print(' ' * (level * 4) + node.name + ':' + str(node.length) + '/ TYPE:' + node.type)
        else:
            print(' ' * (level * 4) + node.support + ':' + str(node.length) + '/ TYPE:' + node.type)
        for child in node.children:
            self.print_tree(child, level + 1)

    def print_node_list(self):
        for node in self.node_list:
            print(node.name)


class RT(Tree):
    tc = None  #Tree Collection
    subtree = None # Selected Subtree

    def __init__(self,name = None,newick_str = None):
        super().__init__(name=name)
        self.parse_newick(newick_str)

    def parse_newick(self,newick_str):
        stack = []
        current_node = None
        name_buffer = ''
        length_buffer = ''
        read_name = False
        read_length = False

        root = Node('root', 0)
        root.type = 'root'
        stack.append(root)

        for char in newick_str:
            if char == '(':
                if self.root is None:
                    self.root = root
                    current_node = Node('')
                    root.children.append(current_node)
                else:
                    new_node = Node('')
                    current_node.children.append(new_node)
                    stack.append(current_node)
                    current_node = new_node
                read_name = True
                read_length = False


            elif char == ')':
                if name_buffer and not current_node.name:
                    if (name_buffer.isdigit()):
                        current_node.name = ''
                        current_node.support = name_buffer
                        current_node.type = INTERNAL
                    else:
                        current_node.name = name_buffer
                        current_node.type = LEAF
                        self.node_list.append(current_node)
                    name_buffer = ''
                if length_buffer:
                    current_node.length = float(length_buffer)
                    length_buffer = ''
                if stack:
                    current_node = stack.pop()
                read_name = True
                read_length = False

            elif char == ',':
                if name_buffer and not current_node.name:
                    if (name_buffer.isdigit()):
                        current_node.name = ''
                        current_node.support = name_buffer
                        current_node.type = INTERNAL
                    else:
                        current_node.name = name_buffer
                        current_node.type = LEAF
                        self.node_list.append(current_node)
                    name_buffer = ''
                if length_buffer:
                    current_node.length = float(length_buffer)
                    length_buffer = ''

                new_node = Node('')
                stack[-1].children.append(new_node)
                current_node = new_node
                read_name = True
                read_length = False

            elif char == ':':
                read_name = False
                read_length = True

            elif read_name:
                name_buffer += char
            elif read_length:
                length_buffer += char

        if name_buffer and not current_node.name:
            if (name_buffer.isdigit()):
                current_node.name = ''
                current_node.support = name_buffer
                current_node.type = INTERNAL
            else:
                current_node.name = name_buffer
                current_node.type = LEAF
                self.node_list.append(current_node)
        if length_buffer:
            current_node.length = float(length_buffer)

        # current_tree.node_list.sort(key=lambda x: x.name, reverse=False)
        # return current_tree

class TC(Tree):
    id = None
    rt = None # Reference Tree's object

    def __init__(self, name = None, id = -1,rt = None):
        super().__init__(name=name)
        self.id = id
        self.rt = rt

        # self.rt.subtree = self


newick_string = "((((((((Ginkgo_biloba:1.00000000000000000000,(Zamia_vazquezii:1.00000000000000000000,(Cycas_micholitzii:1.00000000000000000000,Cycas_rumphii:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(((Cedrus_libani:1.00000000000000000000,Pinus_taeda:1.00000000000000000000)100:1.00000000000000000000,(Prumnopitys_andina:1.00000000000000000000,(Sciadopitys_verticillata:1.00000000000000000000,(Taxus_baccata:1.00000000000000000000,(Juniperus_scopulorum:1.00000000000000000000,Cunninghamia_lanceolata:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(Ephedra_sinica:1.00000000000000000000,(Welwitschia_mirabilis:1.00000000000000000000,Gnetum_montanum:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)98:1.00000000000000000000)100:1.00000000000000000000,(Amborella_trichopoda:1.00000000000000000000,(Nuphar_advena:1.00000000000000000000,(Kadsura_heteroclita:1.00000000000000000000,((Acorus_americanus:1.00000000000000000000,((Yucca_filamentosa:1.00000000000000000000,(Sabal_bermudana:1.00000000000000000000,((Sorghum_bicolor:1.00000000000000000000,Zea_mays:1.00000000000000000000)100:1.00000000000000000000,(Oryza_sativa:1.00000000000000000000,Brachypodium_distachyon:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)96:1.00000000000000000000)29:1.00000000000000000000,(Dioscorea_villosa:1.00000000000000000000,(Colchicum_autumnale:1.00000000000000000000,Smilax_bona-nox:1.00000000000000000000)100:1.00000000000000000000)31:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(((Eschscholzia_californica:1.00000000000000000000,(Aquilegia_formosa:1.00000000000000000000,Podophyllum_peltatum:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(Vitis_vinifera:1.00000000000000000000,(((Boehmeria_nivea:1.00000000000000000000,Medicago_truncatula:1.00000000000000000000)84:1.00000000000000000000,(Larrea_tridentata:1.00000000000000000000,(Populus_trichocarpa:1.00000000000000000000,(Hibiscus_cannabinus:1.00000000000000000000,(Arabidopsis_thaliana:1.00000000000000000000,Carica_papaya:1.00000000000000000000)100:1.00000000000000000000)99:1.00000000000000000000)100:1.00000000000000000000)40:1.00000000000000000000)100:1.00000000000000000000,(Kochia_scoparia:1.00000000000000000000,(Diospyros_malabarica:1.00000000000000000000,((Rosmarinus_officinalis:1.00000000000000000000,(Ipomoea_purpurea:1.00000000000000000000,(Allamanda_cathartica:1.00000000000000000000,Catharanthus_roseus:1.00000000000000000000)100:1.00000000000000000000)98:1.00000000000000000000)100:1.00000000000000000000,(Tanacetum_parthenium:1.00000000000000000000,Inula_helenium:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)72:1.00000000000000000000)80:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(Sarcandra_glabra:1.00000000000000000000,((Persea_americana:1.00000000000000000000,Liriodendron_tulipifera:1.00000000000000000000)100:1.00000000000000000000,(Houttuynia_cordata:1.00000000000000000000,Saruma_henryi:1.00000000000000000000)100:1.00000000000000000000)76:1.00000000000000000000)100:1.00000000000000000000)90:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:2.00000000000000000000,(Equisetum_diffusum:1.00000000000000000000,((Angiopteris_evecta:1.00000000000000000000,(Ophioglossum_petiolatum:1.00000000000000000000,Psilotum_nudum:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(Pteridium_aquilinum:1.00000000000000000000,Alsophila_spinulosa:1.00000000000000000000)100:1.00000000000000000000)70:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,((((Chaetosphaeridium_globosum:1.00000000000000000000,(Coleochaete_irregularis:1.00000000000000000000,Coleochaete_scutata:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(Chara_vulgaris:1.00000000000000000000,(((Uronema_sp:1.00000000000000000000,(Monomastix_opisthostigma:1.00000000000000000000,(Pyramimonas_parkeae:1.00000000000000000000,Nephroselmis_pyriformis:1.00000000000000000000)31:1.00000000000000000000)31:1.00000000000000000000)100:1.00000000000000000000,(Mesostigma_viride:1.00000000000000000000,(Chlorokybus_atmophyticus:1.00000000000000000000,Spirotaenia_minuta:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(Entransia_fimbriata:1.00000000000000000000,Klebsormidium_subtile:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)87:1.00000000000000000000)89:1.00000000000000000000,((Netrium_digitus:1.00000000000000000000,(Roya_obtusa:1.00000000000000000000,(Cosmarium_ochthodes:1.00000000000000000000,Penium_margaritaceum:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(Spirogyra_sp:1.00000000000000000000,(Mesotaenium_endlicherianum:1.00000000000000000000,(Cylindrocystis_brebissonii:1.00000000000000000000,(Cylindrocystis_cushleckae:1.00000000000000000000,Mougeotia_sp:1.00000000000000000000)56:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)6:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,((Nothoceros_aenigmaticus:1.00000000000000000000,Nothoceros_vincentianus:1.00000000000000000000)100:1.00000000000000000000,((Sphagnum_lescurii:1.00000000000000000000,(Polytrichum_commune:1.00000000000000000000,(Physcomitrella_patens:1.00000000000000000000,(Ceratodon_purpureus:1.00000000000000000000,(Hedwigia_ciliata:1.00000000000000000000,((Leucodon_brachypus:1.00000000000000000000,(Thuidium_delicatulum:1.00000000000000000000,(Rhynchostegium_serrulatum:1.00000000000000000000,Anomodon_attenuatus:1.00000000000000000000)93:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(Rosulabryum_cf_capillare:1.00000000000000000000,Bryum_argenteum:1.00000000000000000000)100:1.00000000000000000000)28:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,((Metzgeria_crassipilis:1.00000000000000000000,Bazzania_trilobata:1.00000000000000000000)100:1.00000000000000000000,(Sphaerocarpos_texanus:1.00000000000000000000,(Ricciocarpos_natans:1.00000000000000000000,(Marchantia_emarginata:1.00000000000000000000,Marchantia_polymorpha:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,(Selaginella_moellendorffii_genome:1.00000000000000000000,Selaginella_moellendorffii_1kp:1.00000000000000000000)100:1.00000000000000000000)100:1.00000000000000000000,Huperzia_squarrosa:1.00000000000000000000)100:1.00000000000000000000,Pseudolycopodiella_caroliniana:1.00000000000000000000,Dendrolycopodium_obscurum:1.00000000000000000000);"
# Phylo_Tree = RT(name="RT1",newick_str=newick_string)
# Phylo_Tree = parse_newick(newick_string)
# Phylo_Tree.print_tree()
# Phylo_Tree.print_node_list()
# Phylo_Tree.traverse()