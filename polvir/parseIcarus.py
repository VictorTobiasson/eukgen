# extract summary info from best PDB1-to-PDB2 icarus alignment
def parse_alignment_summary(summary_txtfile):
    import subprocess
    text = subprocess.run(f"cat {summary_txtfile}", shell = True, capture_output=True, text=True).stdout
    # Get first BEST ALIGNMENT section
    section = text.split("BEST ALIGNMENT")[1].split("Aligned distance")[0]
    
    # Initialize collectors for each row type
    collected = {
        'PUs': [],
        'ori. pos': [],
        'connect': [],
        'QUERY': [],
        'match': [],
        'TARGET': [],
        'dist': [],
        'ali. pos': []
    }
    
    current_type = None
    
    for line in section.split('\n'):
        # Skip empty lines
        if not line.strip():
            continue
            
        # Check for row type at the start of line
        for row_type in collected.keys():
            if line.strip().startswith(row_type):
                current_type = row_type
                # Extract the content after the label
                content = line.split(':', 1)[1] if ':' in line else line[len(row_type):].strip()
                collected[current_type].append(content)
                break
        
        # If we're in a current type and line starts with spaces, it's a continuation
        if current_type and line.startswith(' ' * 8) and not any(line.strip().startswith(t) for t in collected.keys()):
            collected[current_type][-1] += line.strip()

    # Join all lines for each type
    result = {key: ''.join(values) for key, values in collected.items()}
    # dictionary, with 'QUERY' and 'dist' keys of interest
    return result

def define_gaps(alignment_data, mindist=4):
    import numpy as np
    window_size = 30
    results = []
    dists = alignment_data["dist"]
    query = alignment_data["QUERY"]
    dists_qinserts = "".join(list(np.array(list(dists))[(np.array(list(query)) != "-") & (np.array(list(query)) != " ")]))
    
    # For each position in string
    for i in range(len(dists_qinserts)):
        pos = dists_qinserts[i]


        isgap_pos = False
        if (pos.isdigit() and int(pos) > mindist) or (pos == " "):
            isgap_pos = True
        # Check all possible windows containing this position
        is_part_of_valid_window = False
        
        # Start checking windows that contain current position
        start_range = max(0, i - window_size + 1)
        end_range = min(i+ 1, len(dists_qinserts) - window_size + 1)
        for start in range(start_range, end_range):
            window = dists_qinserts[start:start + window_size]
            count = sum(1 for c in window if c == ' ' or (c.isdigit() and int(c) > mindist))
            # single possibility of being part of gap = mark as gap
            if isgap_pos:
                if (count / window_size) >= 0.85:
                    is_part_of_valid_window = True
                    break
            # single possibility of not being in a gap = mark as not gap
            elif (count / window_size) < 0.85:
                break
        results.append(is_part_of_valid_window)
    # boolean array, true is gap position
    return results

def map_gaps(bool_arr):
    import numpy as np
    if sum(bool_arr) == 0:
        return [(0,0)] # null gap
    # Add sentinels 
    padded = np.r_[False, bool_arr, False]
    
    # Find transitions
    transitions = np.where(np.diff(padded))[0]
    
    # Get ranges
    return list(zip(transitions[::2], transitions[1::2]))

def crop_pdb(icarus_root, align_root, seq, gap_ranges):
    import os
    from Bio import PDB
    parser = PDB.PDBParser()
    pdb_renum = f"{icarus_root}/{seq}/icarus_results_on_ED_38727.Pavir.Eb01025.1/models/{seq}.renum.pdb"
    structure = parser.get_structure("protein", pdb_renum)
    output_file_tmp = f"{align_root}/{seq}.cropped.tmp.pdb"
    output_file = f"{align_root}/{seq}.cropped.pdb"
    # Create a new structure for the spliced portion
    tmp_structure = PDB.Structure.Structure("spliced_protein")
    # Iterate through models, chains, and residues
    start_number = 5000
    # renumber for each model
    new_resid = start_number
    for model in structure:
        new_model = PDB.Model.Model(model.id)
        tmp_structure.add(new_model)
        for chain in model:
            new_chain = PDB.Chain.Chain(chain.id)
            new_model.add(new_chain)
            # new_resid = start_number
            for residue in chain:
                # check if insertion or not
                residue_num = residue.id[1]
                if not any(start <= residue_num - 1 < end for start, end in gap_ranges):
                    # renum
                    residue.id = (residue.id[0], new_resid, residue.id[2])
                    new_resid += 1
                    new_chain.add(residue)
    io = PDB.PDBIO()
    io.set_structure(tmp_structure)
    io.save(output_file_tmp)
    
    # Create a new structure for the spliced portion
    structure = parser.get_structure("protein", output_file_tmp)
    final_structure = PDB.Structure.Structure("spliced_protein")
    # Iterate through models, chains, and residues
    start_number = 1
    # renumber for each model
    new_resid = start_number
    for model in structure:
        new_model = PDB.Model.Model(model.id)
        final_structure.add(new_model)
        for chain in model:
            new_chain = PDB.Chain.Chain(chain.id)
            new_model.add(new_chain)
            # new_resid = start_number
            for residue in chain:
                residue.id = (residue.id[0], new_resid, residue.id[2])
                new_chain.add(residue)
                new_resid += 1
    io = PDB.PDBIO()
    io.set_structure(final_structure)
    io.save(output_file)
    os.system(f"rm {output_file_tmp}")
    return output_file
def entropy_filter(basepath, threshold=0.2):
    import os
    from core_functions.helper_functions import fasta_to_dict, dict_to_fasta, filter_by_entropy 
    for f in os.listdir(basepath):
        if f.endswith((".fasta", ".fa", ".afa")):
            inpath = f"{basepath}{f}"
            file_ext = f.split(".")[-1]
            f_out = f.replace(file_ext, ".entropy_filt.fasta")
            outpath = f"{basepath}{f_out}"
            poly_dict = fasta_to_dict(file=inpath)
            poly_filter = filter_by_entropy(poly_dict, threshold)
            dict_to_fasta(poly_filter, write_file=outpath)