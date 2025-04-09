version 1.0

workflow NM_Analysis_Completion {
    input {
        File assembly
        String sample_id
    }

    call Get_Serogroup {
        input: 
            assembly=assembly,
            sample_id=sample_id

    }

    call Get_Clonal_Complex {
        input: assembly=assembly
    }

    call Get_Mendevar {
        input: assembly=assembly
    }

    output {
        String clonal_complex = Get_Clonal_Complex.clonal_complex
        String mendevar_bexsero_reactivity = Get_Mendevar.bexsero
        String mendevar_trumenba_reactivity = Get_Mendevar.trumenba
        String serogroup = Get_Serogroup.serogroup
        String genogroup = Get_Serogroup.genogroup
    }
}

task Get_Serogroup {
    input {
        String sample_id
        File assembly
    }

    command <<<
        python3 /app/nm_completion.py \
        --input ~{assembly} \
        --serogroup_filename ~{sample_id}_serogroup.txt \
        --genogroup_filename ~{sample_id}_genogroup.txt \
        --sample ~{sample_id}
    >>>

    output {
        String serogroup = read_string("~{sample_id}_serogroup.txt")
        String genogroup = read_string("~{sample_id}_genogroup.txt")
    }

    runtime {
        docker: "bioinfomoh/nm_completion_analysis:1"
    }
}

task Get_Clonal_Complex {
    input {
        File assembly
    }

    command <<<
        url="https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/1/sequence"

        (echo -n '{"base64":true,"sequence": "'; base64 ~{assembly}; echo '"}') | \
        curl -s -H "Content-Type: application/json" -X POST ${url} -d @- | \
        jq -r '.fields.clonal_complex' > cc.txt
    >>>

    output {
        String clonal_complex = read_string("cc.txt")
    }

    runtime {
        docker: "devorbitus/ubuntu-bash-jq-curl:latest"
    }
}

task Get_Mendevar {
    input {
        File assembly
    }

    command <<<
        URL_BAST="https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/53/sequence"
        URL_BAST_DESIGNATION="https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/53/designations"

        # Send the sequence to the server and save the response to bast.json
        (echo -n '{"base64":true,"sequence": "'; base64 ~{assembly}; echo '"}') | \
        curl -s -H "Content-Type: application/json" -X POST ${URL_BAST} -d @- | \
        jq . > bast.json

        declare -A loci
        # Read the content of bast.json and extract allele_id for each key in 'exact_matches'
        for key in $(jq -r '.exact_matches | keys_unsorted | .[]' bast.json); do
            allele_id=$(jq -r ".exact_matches[\"$key\"][0].allele_id" bast.json)
            loci["${key}"]=$allele_id
        done

        data=("fHbp_peptide" "PorA_VR2" "PorA_VR1" "NHBA_peptide" "NadA_peptide")

        url=""
        for val in "${data[@]}"; do
            if [[ -n "${loci[$val]}" ]]; then  # if the key exists in "al"
                url+='"'$val'":[{"allele":"'"${loci[$val]}"'"}],'

            else
                url+='"'$val'":[{"allele":"'"0"'"}],'    
            fi
        done

        url="${url%?}"

        curl -s -H "Content-Type: application/json" \
        -X POST "${URL_BAST_DESIGNATION}" \
        -d "{\"designations\": { ${url} }}" > bast_type.json


        # Extract Bexsero and Trumenba reactivity
        bexsero=$(jq -r '.fields.MenDeVAR_Bexsero_reactivity' bast_type.json > bexsero.txt)
        trumenba=$(jq -r '.fields.MenDeVAR_Trumenba_reactivity' bast_type.json > trumenba.txt)

    >>>

    output {
        File bast_type = 'bast_type.json'
        String bexsero = read_string("bexsero.txt")
        String trumenba = read_string("trumenba.txt")
    }

    runtime {
        docker: "devorbitus/ubuntu-bash-jq-curl:latest"
    }
}

