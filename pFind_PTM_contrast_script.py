# coding = utf-8
import copy, os

wd_path = r"L:\pFind\pFind_result\ID_data_Hela_QE_HF_120min\Comparison" # important the pFind.protein path
output_name = "pFind_PTM_contrast_result.txt" # important! the output file name
modify = "Oxidation[M]" # target modification name


def get_modi_info(linelist, mod_tgt):
    pros = linelist[10]
    mods = linelist[8]
    pep = linelist[3]
    pros_list = pros[:-1].split('/')
    pep_pos_list = linelist[11][:-1].split('/')
    mod_list = mods[:-1].split(";")
    spec_num = linelist[-1]
    score = linelist[7]
    # spec_title = linelist[-3]
    info_list = []
    for mod_info in mod_list:
        md_pos, md_name = mod_info.split(',')
        if md_name  == mod_tgt:
            pro_pos_list = []
            for i in range(len(pros_list)):
                pro = pros_list[i]
                pep_pos = int(pep_pos_list[i].split(',')[0])
                mod_pos = pep_pos + int(md_pos)
                pro_pos = pro + "[" + str(mod_pos) + "]"
                pro_pos_list.append(pro_pos)
            pro_poss = "/".join(pro_pos_list) + "/"
            info_list.append([pros, pro_poss, pep, mods, score, spec_num])
    
    return info_list

                
def get_complement_proteins(protein_list, i):
    copy_list = copy.deepcopy(protein_list)
    del copy_list[i]
    if copy_list == []:
        return ""
    else:
        return "/".join(copy_list)


def get_info_dic(fl_path):
    uni_dic = {}

    f = open(fl_path).readlines()
    for line in f[:]:
        if modify in line:
            linelist = line.rstrip("\n").split("\t")
            # print(linelist)
            pro_sites_info = get_modi_info(linelist, modify)
            # print()
            for pro_site in pro_sites_info:
                pros, pro_poss, pep, mods, score, spec_num = pro_site
                
                protein_list = pros[:-1].split("/")
                is_uni = "unique"
                if len(protein_list) > 1:
                    is_uni = "not unique"
                # pep_info_list.append(is_uni)
                pro_poss_list = pro_poss[:-1].split('/')

                for i in range(len(protein_list)):
                    pep_info_list = [pep, mods, score, spec_num]
                    pep_info_list.append(is_uni)
                    protein = protein_list[i]
                    complem_proteins = get_complement_proteins(protein_list, i)
                    pep_info_list.append(complem_proteins)
                    protein_pos = pro_poss_list[i]
                    if protein not in uni_dic:
                        uni_dic[protein] = {}
                        uni_dic[protein][protein_pos] = [pep_info_list]
                    else:
                        if protein_pos not in uni_dic[protein]:
                            uni_dic[protein][protein_pos] = [pep_info_list]
                        else:
                            uni_dic[protein][protein_pos].append(pep_info_list)
    return uni_dic


def reformate_dic(uni_dic):
    # pros_sites = []
    # for pros in uni_dic:
    #     pros_sites.append([pros, len(uni_dic[pros])])
    
    # new_pros_seq = [x[0] for x in sorted(pros_sites, key = lambda x:x[1], reverse= True)]
    re_uni_dic = {}
    for pros in uni_dic:
        pro_m_num = 0
        pro_best_score = 1
        pro_p_num = 0
        pro_site_num = 0
        poss_pep_dic = uni_dic[pros]
        # wline_list = []
        have_uni_pep_pro = "no unique peptides"
        pro_sites_dic = {}
        for site in poss_pep_dic:
            w_sub_line = []
            site_m_num = 0
            site_best_score = 1
            site_unique_best_score = 1
            site_p_num = 0
            unique_pep = 0
            unique_spec= 0
            have_unipue_pep_site = "no unique peptide"
            site_peps = poss_pep_dic[site]
            pep_info_dic = {}
            for pep_info in site_peps:
                pep, mods, score, spec_num, is_uni, shared_pros = pep_info
                pep_info_dic[(pep, mods)] = [spec_num, shared_pros, score, is_uni]
                site_m_num += int(spec_num)
                site_p_num += 1
                if float(score) < site_best_score:
                    site_best_score = float(score)
                if is_uni == "unique":
                    have_unipue_pep_site = "have unique peptides"
                    have_uni_pep_pro = "have unique peptides"
                    unique_pep += 1
                    unique_spec += int(spec_num)
                    if float(score) < site_unique_best_score:
                        site_unique_best_score = float(score)
            pro_sites_dic[site] = [[str(site_p_num) + "(" + str(unique_pep) + ")", str(site_m_num) + "(" + str(unique_spec) + ")", site_best_score, have_unipue_pep_site], pep_info_dic]
         
            pro_m_num += site_m_num
            pro_p_num += site_p_num
            pro_site_num += 1
            if pro_best_score > site_best_score:
                pro_best_score = site_best_score
        re_uni_dic[pros] = [[pro_site_num, pro_p_num, pro_m_num, have_uni_pep_pro], pro_sites_dic]
        
    return re_uni_dic


def get_all_pros_sites_peps(total_dic):
    struced_dic = {}
    for sp in total_dic:
        re_uni_dic = total_dic[sp]
        for pro in re_uni_dic:
            if pro not in struced_dic:
                struced_dic[pro] = [[], {}]
            pro_sites_dic = re_uni_dic[pro][1]
            
            for site in pro_sites_dic:
                if site not in struced_dic[pro][1]:
                    struced_dic[pro][1][site] = [[], {}]
        
                site_peps_dic = pro_sites_dic[site][1]
                
                for pep_m in site_peps_dic:
                    struced_dic[pro][1][site][1][pep_m] = []
    return struced_dic


def compare_results(total_dic):
    if len(total_dic) == 1:
        return list(total_dic.values())[0]
    else:
        structed_dic = get_all_pros_sites_peps(total_dic)
        sp_list = list(total_dic.keys())
        for pro in structed_dic:
            for sp in sp_list:
                uni_dic = total_dic[sp]
                if pro not in uni_dic:
                    structed_dic[pro][0].extend(["", "", "", ""])
                    for site in structed_dic[pro][1]:
                        structed_dic[pro][1][site][0].extend(["", "", "", ""])

                        for pep_m in structed_dic[pro][1][site][1]:
                            structed_dic[pro][1][site][1][pep_m].extend(["", "", "", ""])
                else:
                    structed_dic[pro][0].extend(uni_dic[pro][0])
                    for site in structed_dic[pro][1]:
                        if site not in uni_dic[pro][1]:
                            structed_dic[pro][1][site][0].extend(["", "", "", ""])
                            
                            for pep_m in structed_dic[pro][1][site][1]:
                                structed_dic[pro][1][site][1][pep_m].extend(["", "", "", ""])

                        else:
                            structed_dic[pro][1][site][0].extend(uni_dic[pro][1][site][0])
                            for pep_m in structed_dic[pro][1][site][1]:
                                if pep_m not in uni_dic[pro][1][site][1]:
                                    structed_dic[pro][1][site][1][pep_m].extend(["", "", "", ""])
                                else:
                                    structed_dic[pro][1][site][1][pep_m].extend(uni_dic[pro][1][site][1][pep_m])
        return structed_dic


def write_dic2_file(final_dic, sp_list, w_name):
    b = open(os.path.join(wd_path, w_name), 'w')
    line1_list = [""] * 4
    for sp in sp_list:
        line1_list.extend([sp, "", "", ""])
    b.write("\t".join(line1_list)+"\n")
    line2_list = ["Protein", "", "", ""]
    line2_list.extend(["Total_site_num@pro","Total_pep_num@pro", "Total_spec_num@pro", "have unique peptides?"]* len(sp_list))
    b.write("\t".join(line2_list)+"\n")
    line3_list = ["", "site", "", ""]
    line3_list.extend(["Total_pep_num@site(Total_unique_pep_num@site)", "Total_spec_num@site(Total_unique_spec_num@site)", "best-score@site", "have unique peptides?"]* len(sp_list))
    b.write("\t".join(line3_list)+"\n")
    line4_list = ["", "", "Peptide", "Modification"]
    line4_list.extend(["Total_spec_num@pep", "shared proteins", "best-score@pep", "is unique?"]* len(sp_list))
    b.write("\t".join(line4_list)+"\n")

    for pro in final_dic:
        pro_info, sites_dic = final_dic[pro]
        pro_wlist = [pro, "", "", ""]
        pro_wlist.extend(pro_info)
        b.write("\t".join([str(x) for x in pro_wlist]) + "\n")
        for site in sites_dic:
            site_info, peps_dic =  sites_dic[site]
            site_wlist = ["", site, "", ""]
            site_wlist.extend(site_info)
            b.write("\t".join([str(x) for x in site_wlist]) + "\n")
            for pep, mod in peps_dic:
                pep_wlist = ["", "", pep, mod]
                pep_wlist.extend(peps_dic[pep,mod])
                b.write("\t".join([str(x) for x in pep_wlist]) + "\n")
    b.close()



def main():
    total_dic = {}
    for fl in os.listdir(wd_path):
        if fl.endswith(".protein"):
            fl_path = os.path.join(wd_path, fl)
            fl_name = fl[:-8]
            total_dic[fl_name] =  reformate_dic(get_info_dic(fl_path))
    final_dic = compare_results(total_dic)
    sp_list = list(total_dic.keys())
    print(final_dic)
    write_dic2_file(final_dic, sp_list, output_name)
    


if __name__ == "__main__":
    main()
    print("Well Done!")