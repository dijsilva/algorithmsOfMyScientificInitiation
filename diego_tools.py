import matplotlib.pyplot as plt

def blast_out_Normal(path_ncbi_info, path_pre_tax):
    ncbi = open(path_ncbi_info, 'r')
    tax = open(path_pre_tax, 'r')
    ncbi_info = []
    pre_tax = []
    
    #as informações vindas destes banco de dados não obedecem bem um padrão, portanto foram retirados
    ignore = ['pdb', 'sp', 'prf']
    for info in ncbi:
        ncbi_info.append(info.split("\n")[0 ])
    for info2 in tax:
        pre_tax.append(info2.split("\n")[0])

    #sel_ncbi é uma lista que conterá os códigos das sequencias que apresentarem as caracteristicas compativeis
    sel_ncbi = []
    for pos, ptn in enumerate(pre_tax):
        #checar a formatação do output do blast
        if len(ncbi_info[pos].split(' ')[0].split('|')) > 4:
            #checar se o banco de dados da sequencia não está na lista de banco de dados problemáticos e há os colchetes no arquivo de taxonomia (é indicativo de há o nome de uma espécia dentro do output)
            if ncbi_info[pos].split(' ')[0].split('|')[2] not in ignore:
                if '[' in ptn:
                    if ']' in ptn:
                        sel_ncbi.append(pos)

    org = []
    for posi in sel_ncbi:
        #splitar as especies do output do blast para os hits selecionados na etapa anterior deste algoritmo
        org_temp = []
        selection = pre_tax[posi]
        selection = selection.split('>')
        
        if len(selection) > 1:
            #percorrer a lista com as espécies
            for sel in selection:
                if '[' in sel:
                    if ']' in sel:
                        init = sel.find('[')
                        end = sel.find(']')
                        #retirar o nome da espécie
                        organism = sel[(init+1): end]
                        organism = organism.split(' ')
                        #salvando a espécie na lista temporaria
                        if len(organism) > 1:
                            organism = ("{} {}".format(organism[0], organism[1]))
                            org_temp.append(organism)
                        else:
                            organism = ("{}".format(organism[0]))
                            org_temp.append(organism)
        else:
            sel = selection[0]
            if '[' in sel:
                if ']' in sel:
                    init = sel.find('[')
                    end = sel.find(']')
                    organism2 = sel[(init+1): end]
                    organism2 = organism2.split(' ')
                    if len(organism2) > 1:
                        organism2 = ("{} {}".format(organism2[0], organism2[1]))
                        org_temp.append(organism2)
                    else:
                        organism2 = ("{}".format(organism2[0]))
                        org_temp.append(organism2)
        
        #ordenando a lista para retirar as espécie repetidas
        org_temp = set(org_temp)
        org_temp = sorted(list(org_temp))
        org.append(org_temp)

    species = {}
    for lists in org:
        for sp in lists:
            if sp not in species:
                species[sp] = 1
            else:
                species[sp] = species[sp] + 1

    sorted_species = sorted(species, key=species.get)
    tuple_of_species  = sorted(species.items(), key=lambda x:x[1])

    count_species = {}
    for lists in tuple_of_species:
        if lists[1] not in count_species:
            count_species[lists[1]] = 1
        else:
            count_species[lists[1]] = count_species[lists[1]] + 1
    return ncbi_info, pre_tax, sel_ncbi, org, count_species, sorted_species, tuple_of_species


