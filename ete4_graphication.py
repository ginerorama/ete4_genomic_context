import sys
from ete4 import Tree, random_color
from ete4.smartview import TreeLayout, ArrowFace, TextFace
from ete4.smartview.renderer.layouts.default_layouts import LayoutLeafName
from ete4.smartview.renderer.layouts.context_layouts import LayoutGenomicContext


def get_layout_context_old(width=150, height=15):
    def layout_fn(node):
        if node.is_leaf():
            context = node.props.get("_context")
        else:
            first_leaf = next(node.iter_leaves())
            context = first_leaf.props.get("_context")
        if context:
            for idx, gene in enumerate(context):
                name = gene.get("name")
                color = gene.get("color", "gray")
                strand = gene.get("strand", "+")
                orientation = "left" if strand == "-1" else "right"
                cluster = gene.get("cluster")
                description = gene.get("description")
                domains = gene.get("domains")
                orf = gene.get("orf")
                genome_name = gene.get("genome_name")
                text = name 
                props = {"ORF": orf, "Gene_symbol": name, "Cluster": cluster, "Description": description,"Prosite domains":domains ,"Genome":genome_name}
                #tooltip = ["Gene information\n"]
                link = "<p>Crea un enlace a\n<a href='https://www.mozilla.org/es-ES/'>la página de inicio de Mozilla</a>.\n</p>\n"
                tooltip = [link]
                tooltip.append("\n".join([ key+": "+value for key,value in props.items()]))
                print(tooltip)
                if idx == nside: 
                    stroke_color = "red"
                    stroke_width = "1.5px"
                else:
                    stroke_color = "gray"
                    stroke_width = "1px"
         
                arrow = ArrowFace(width=width, height=height,
                        orientation=orientation, color=color,
                        text=text, tooltip=tooltip,max_fsize=11,stroke_color=stroke_color, stroke_width=stroke_width,
                        padding_x=1.5, padding_y=1.5)
                text = TextFace(name, padding_x=2, padding_y=2,max_fsize = 6)
                node.add_face(arrow, position="aligned", column=idx,
                        collapsed_only=(not node.is_leaf()))
                # node.add_face(text, position="aligned", column=idx,
                        # collapsed_only=(not node.is_leaf()))
   
    return layout_fn



def set_tree_style_old(style):
    style.collapse_size =1


def get_context(path, nside,eggnog_mapper, file):
    data = []
    if eggnog_mapper == True:
        eggnogmapper_hash = read_eggnogmapper_file(file)
    with open("./gcontext.csv") as handle:
        headers = handle.readline().strip().split(",")
        for line in handle.readlines():
            orf, cluster, strand, description, gname, genome_name,domains,sequences,conservation_score = line.strip().split(",")
            domains = "\n"+domains.replace("@","\t").replace("|","\t")
            if eggnog_mapper == True:
                if orf in eggnogmapper_hash:
                    eggmapper= [v for k,v in eggnogmapper_hash[orf].items() ]
                    eggNOG_OGs,COG_category,eggnog_description,Preferred_name,KEGG_ko,PFAMs = eggmapper
                    gene = { "orf": orf, "cluster": cluster,"eggNOG_OGs":eggNOG_OGs.replace(","," "),"COG_category":COG_category,
                     "strand": strand, "Eggnog_description":eggnog_description,"KEGG_ko":KEGG_ko.replace("ko:",""),"PFAMs":PFAMs,
                     "name": gname if gname != ""else Preferred_name.split(",")[0], "description": description,  
                     "genome_name":genome_name, "domains":domains, 
                     "conservation_score": conservation_score}
                else:
                    gene = { "orf": orf, "cluster": cluster,"strand": strand, "name": gname, "description": description,  
                     "genome_name":genome_name, "domains":domains,"conservation_score": conservation_score}
            
            else:
                gene = { "orf": orf, "cluster": cluster,"strand": strand, "name": gname, "description": description,  
                     "genome_name":genome_name, "domains":domains,"conservation_score": conservation_score}
      
            if gene["name"] == " - ":
                gene["name"] == ""
            data.append(gene)
            #if orf == "FOC17_RS22470":
                #print(gene)

    window = 2 * nside + 1
    context = {}
    genome_ID_for_anchor = {}
    for idx in range(1, len(data) // window + 1):
        anchor = idx * window - nside - 1
        context[data[anchor]["orf"]] = data[anchor-nside:anchor+nside+1]
        genome_ID = data[anchor]["genome_name"]
        genome_ID_for_anchor[data[anchor]["orf"]] = genome_ID

    #print(context)
    return context,genome_ID_for_anchor


# Get colors from unique "clusters"
def get_color_from_unique_clusters(context,conservation_threshols):
    unique_clusters = set()
    for v in context.values():
        #unique_clusters.update(gene["cluster"] for g in v)
        for gene in v:
            if float(gene["conservation_score"]) > float(conservation_threshols):
                #print(gene["conservation_score"],float(conservation_threshols))
                unique_clusters.add(gene["cluster"])
            else:
                gene["cluster"] = '0'
                unique_clusters.add(gene["cluster"])

    if '0' in unique_clusters:
        unique_clusters.remove('0')
        colors = random_color(num=len(unique_clusters), l=0.5, s=0.5)
        cluster2color = dict(zip(unique_clusters, colors))
        cluster2color['0']= "#d0d0d0" #If pseudogen we add grey color, change the face in future by broken arrow?
    else:
        colors = random_color(num=len(unique_clusters), l=0.5, s=0.5)
        cluster2color = dict(zip(unique_clusters, colors))

    for v in context.values():
        for g in v:
            g["color"] = cluster2color[g["cluster"]]


def read_eggnogmapper_file(file):
	"""
	Eggnogmapper fields in csv:

	query
	seed_ortholog
	evalue
	score
	eggNOG_OGs
	max_annot_lvl
	COG_category
	Description
	Preferred_name
	GOs
	EC
	KEGG_ko
	KEGG_Pathway
	KEGG_Module
	KEGG_Reaction
	KEGG_rclass
	BRITE
	KEGG_TC
	CAZy
	BiGG_Reaction
	PFAMs

	"""

	eggmapper = open(file,"r")

	eggmapper_hash =  {}

	for line in eggmapper:
		if not line.startswith("#"):

			[query,
			seed_ortholog,
			evalue,
			score,
			eggNOG_OGs,
			max_annot_lvl,
			COG_category,
			Description,
			Preferred_name,
			GOs,
			EC,
			KEGG_ko,
			KEGG_Pathway,
			KEGG_Module,
			KEGG_Reaction,
			KEGG_rclass,
			BRITE,
			KEGG_TC,
			CAZy,
			BiGG_Reaction,
			PFAMs] = line.rstrip().split("\t")

			eggmapper_hash[query] = {"eggNOG_OGs":eggNOG_OGs,"COG_category":COG_category,"Description":Description,
			"Preferred_name":Preferred_name,"KEGG_ko":KEGG_ko,"PFAMs":PFAMs}

	return eggmapper_hash


nside = 3   

if sys.argv[2]:
    eggnog_mapper = True
    file = sys.argv[2]
    context, genome_ID_for_anchor = get_context("./gcontext.csv", nside, eggnog_mapper,file)
else:
    eggnog_mapper = False
    file = ""
    context, genome_ID_for_anchor = get_context("./gcontext.csv", nside, eggnog_mapper,file)
    pass

#context, genome_ID_for_anchor = get_context("./gcontext.csv", nside, eggnog_mapper)
get_color_from_unique_clusters(context, 0)    
#print(context)

tree_file = open(sys.argv[1],"r").readline()
t = Tree(tree_file)
#print(t)
#t.populate(len(anchors))

for leaf in t.iter_leaves():
    spname = genome_ID_for_anchor[leaf.name]
    leaf.props["_context"] = context[leaf.name]  # add context data to leaf (anchor)
    leaf.props["gene"] = leaf.name 
    leaf.name += " ("+spname+")"
    

leafname = LayoutLeafName(max_fsize=10)
leafname.active = True


#for n in t.traverse():
#    n.dist = 0.5

layouts = []
for v in [0,10,25,50,75]:
    context_layout = LayoutGenomicContext("conservation_"+str(v),conservation_threshold = v, nside = nside, collapse_size = 1)
    context_layout.active = v == 0
    layouts.append(context_layout)

layouts.append(leafname)
t.explore("context", layouts = layouts)


