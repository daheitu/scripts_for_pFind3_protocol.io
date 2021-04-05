# scripts_for_pFind3_protocol.io
The scripts are used for compare two or more identification results from **pFind 3**. please follow the steps:
1. Create a new folder in a desired location and give it a name. 
2. To compare across samples either the identified proteins or PTMs, respectively, copy either “pFind_protein_contrast_script.py” or “pFind_PTM_contrast_script.py” and the different pFind.protein files to be compared into this new folder. Note that each pFind.protein file should be renamed before it comes to this folder. The file name should be a ready reminder of the sample and the purpose of the search. 
3. Open “pFind_protein_contrast_script.py” or “pFind_PTM_contrast_script.py” using ‘Notepad++’ or another editor, specify the path of the new folder and save. You may change the name of the output file if you dislike the default one. For pFind_PTM_contrast_script.py, one particular modification must be specified.
4. Open ‘Windows Explorer’, type ‘cmd’ in the location bar, and press ‘Enter’ to open ‘Command Prompt’. 
5. Lastly, run the script by entering the following on the command line:<br>
`python pFind_protein_contrast_script.py`<br>
Or, for PTMs:<br>
`python pFind_PTM_contrast_script.py`
6. Open with Excel the contrast file generated from the last step and explore. Header annotation of pFind_protein_contrast_result.txt and that of or pFind_PTM_contrast_result.txt can be found in Table S9 and Table S10 in SI. 
