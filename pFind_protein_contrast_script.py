# coding = utf-8
import copy, os
# from cyMStools import plink2_summary

wd_path = r"L:\pFind\pFind_result\ID_data_Hela_QE_HF_120min\Comparison" # important the pFind.protein path
output_name = "pFind_protein_contrast_result.txt" # important! the output file name


#求互补蛋白序列                
def get_complement_proteins(protein_list, pro, group_pro):
    if len(protein_list) == 1 and protein_list[0] == pro:
        return ""
    else:
        group_set = set(group_pro) # group proteins
        pep_related_pro = set(protein_list)
        target_pro_set = set(pro)
        complement_pros = pep_related_pro - target_pro_set & group_set
        return "/".join(list(complement_pros))

    
def get_info_dic(fl_path):
    rep_dic = {}

    f = open(fl_path).readlines()
    i = 0
    while i < len(f)-5:
        if f[i].startswith("-----"):
            break
        else:
            linelist = f[i].rstrip("\n").split("\t")
            if not linelist[0].isdigit():
                i += 1
            else:
                print(i)
                group_pro = [linelist[1]]
                i += 1
                while i < len(f)-5:
                    sub_list = f[i].rstrip("\n").split("\t")
                    if sub_list[2].isdigit():
                        break
                    else:
                        group_pro.append(sub_list[2])
                        i += 1
                while i < len(f)-5:
                    pep_info_list = f[i].rstrip("\n").split("\t")
                    if "-----" in f[i]:
                        i = len(f)
                    elif pep_info_list[0].isdigit():
                        break
                    else:
                        pep = pep_info_list[3]
                        mod = pep_info_list[8]
                        pros = pep_info_list[10]
                        spec = pep_info_list[-1]
                        pro_list = pros[:-1].split("/")
                        best_score = pep_info_list[7]
                        is_unique = "not unique"
                        if len(pro_list) == 1:
                            is_unique = "is unique"
                        for pro in pro_list:
                            if pro in group_pro:
                                comple_pros = get_complement_proteins(pro_list, pro, group_pro)
                                if pro not in rep_dic:
                                    rep_dic[pro] = {}
                                    rep_dic[pro][pep, mod] = [spec, best_score, is_unique, comple_pros]
                                else:
                                    rep_dic[pro][pep, mod] = [spec, best_score, is_unique, comple_pros]

                        i += 1
    return rep_dic


def reformate_dic(uni_dic):
    re_uni_dic = {}
    for pro in uni_dic:
        pro_p_num = 0
        pro_pep_unique = 0
        pro_spec_num = 0
        pro_spec_unique = 0
        poss_pep_dic = uni_dic[pro]
        have_uni_pep_pro = "no unique peptide"
        for pep_info in poss_pep_dic:
            spec, best_score, is_unique, comple_pros = poss_pep_dic[pep_info]
            pro_p_num += 1
            pro_spec_num += int(spec)
            if is_unique == "is unique":
                pro_pep_unique += 1
                pro_spec_unique += int(spec)
                have_uni_pep_pro= "have unique peptide"
        re_uni_dic[pro] = [["%d(%d)" % (pro_p_num, pro_pep_unique), "%d(%d)" % (pro_spec_num, pro_spec_unique),have_uni_pep_pro], poss_pep_dic]
        
    return re_uni_dic


def get_all_pros_sites_peps(total_dic):
    struced_dic = {}
    for sp in total_dic:
        re_uni_dic = total_dic[sp]
        for pro in re_uni_dic:
            if pro not in struced_dic:
                struced_dic[pro] = [[], {}]
            pro_pep_dic = re_uni_dic[pro][1]
            
            for pep in pro_pep_dic:
                if pep not in struced_dic[pro][1]:
                    struced_dic[pro][1][pep] = []
        
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
                    for pep in structed_dic[pro][1]:
                        structed_dic[pro][1][pep].extend(["", "", "", ""])

                else:
                    structed_dic[pro][0].extend(uni_dic[pro][0])
                    structed_dic[pro][0].append("")
                    for pep in structed_dic[pro][1]:
                        if pep not in uni_dic[pro][1]:
                            structed_dic[pro][1][pep].extend(["", "", "", ""])
                            
                        else:
                            structed_dic[pro][1][pep].extend(uni_dic[pro][1][pep])
        return structed_dic


def write_dic2_file(final_dic, sp_list, w_name):
    b = open(os.path.join(wd_path, w_name), 'w')
    line1_list = [""] * 3
    for sp in sp_list:
        line1_list.extend([sp, "", "", ""])
    b.write("\t".join(line1_list)+"\n")
    line2_list = ["Protein", "", ""]
    line2_list.extend(["Total_pep_num@pro(Total_unique_pep_num@pro)","Total_spec_num@pro(Total_unique_spec_num@pro)", "have unique peptide?", ""]* len(sp_list))
    b.write("\t".join(line2_list)+"\n")
    # line3_list = ["", "site", "", ""]
    # line3_list.extend(["Total_pep_num@site(Total_unique_pep_num@site)", "Total_spec_num@site(Total_unique_spec_num@site)", "best-score@site", "have unique peptide?"]* len(sp_list))
    # b.write("\t".join(line3_list)+"\n")
    line4_list = ["", "Peptide", "Modification"]
    line4_list.extend(["Total_spec_num@pep", "best-score@pep", "is unique?", "shared proteins"]* len(sp_list))
    b.write("\t".join(line4_list)+"\n")

    for pro in final_dic:
        pro_info, pep_dic = final_dic[pro]
        pro_wlist = [pro, "", ""]
        pro_wlist.extend(pro_info)
        b.write("\t".join([str(x) for x in pro_wlist]) + "\n")
        for pep, mod in pep_dic:
            pep_wlist = ["", pep, mod]
            pep_wlist.extend(pep_dic[pep,mod])
            b.write("\t".join([str(x) for x in pep_wlist]) + "\n")
    b.close()



def main():
    total_dic = {}
    for fl in os.listdir(wd_path):
        if fl.endswith(".protein"):
            fl_path = os.path.join(wd_path, fl)
            fl_name = fl[:-8]
            print(fl_path)
            total_dic[fl_name] =  reformate_dic(get_info_dic(fl_path))
    final_dic = compare_results(total_dic)
    sp_list = list(total_dic.keys())
    # print(final_dic)
    write_dic2_file(final_dic, sp_list, output_name)
    


if __name__ == "__main__":
    main()
    print("Well Done!")
