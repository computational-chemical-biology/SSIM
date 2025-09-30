import pyopenms as oms
import os

def getFeatureMaps(mzML_list):
    feature_maps = []
    for file in mzML_list:
        # carregar os arquivos mzML file no objeto MSExperiment
        print('Processing file:', file) 
        exp = oms.MSExperiment()
        oms.MzMLFile().load(
            file, exp
        )  
    
        # detecção de 'traços'  de massa
        mass_traces = (
            []
        )  # criar lista vazia onde as massas serão armazenadas
        mtd = oms.MassTraceDetection()
        mtd_par = (
            mtd.getDefaults()
        )  # obter parâmetros padrão para edição
        mtd_par.setValue("mass_error_ppm", 15.0)  # instrumentos de alta resolução, exemplo: orbitraps
        mtd_par.setValue(
            "noise_threshold_int", 1.0e03
        )  # dependente dos dados (usualmente funciona para orbitraps)
        mtd.setParameters(mtd_par)  # definir os novos parâmetros
        mtd.run(exp, mass_traces, 0)  # executar detecção de massas
    
        # detecção de picos de eluição
        mass_traces_deconvol = []
        epd = oms.ElutionPeakDetection()
        epd_par = epd.getDefaults()
        epd_par.setValue(
            "width_filtering", "fixed"
        )  # Os filtros de parâmetros fixos removem traços fora do intervalo [min_fwhm: 1.0, max_fwhm: 60.0]
        epd.setParameters(epd_par)
        epd.detectPeaks(mass_traces, mass_traces_deconvol)
    
        # Detecção de features
        feature_map = oms.FeatureMap()  # features encontradas
        chrom_out = []  # lista de cromatogramas
        ffm = oms.FeatureFindingMetabo()
        ffm_par = ffm.getDefaults()
        ffm_par.setValue(
            "remove_single_traces", "true"
        )  # remove traços sem traços isotópicos satélite
        ffm.setParameters(ffm_par)
        ffm.run(mass_traces_deconvol, feature_map, chrom_out)
        feature_map.setUniqueIds()  # Atribui uma nova, id única válida por feature
        feature_map.setPrimaryMSRunPath(
            [file.encode()]
        )  # Ajusta o caminho do arquivo para o local dos aquivos mzML
        feature_maps.append(feature_map)
    
    return feature_maps

def alignFeatureRT(feature_maps):
    # use as reference for alignment, the file with the largest number of features
    # (works well if you have a pooled QC for example)
    ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])
    
    aligner = oms.MapAlignmentAlgorithmPoseClustering()
    
    trafos = {}
    
    # parameter optimization
    aligner_par = aligner.getDefaults()
    aligner_par.setValue("max_num_peaks_considered", -1)  # infinite
    aligner_par.setValue(
        "pairfinder:distance_MZ:max_difference", 10.0
    )  # Never pair features with larger m/z distance
    aligner_par.setValue("pairfinder:distance_MZ:unit", "ppm")
    aligner.setParameters(aligner_par)
    aligner.setReference(feature_maps[ref_index])
    
    for feature_map in feature_maps[:ref_index] + feature_maps[ref_index + 1 :]:
        trafo = oms.TransformationDescription()  # save the transformed data points
        aligner.align(feature_map, trafo)
        trafos[feature_map.getMetaValue("spectra_data")[0].decode()] = trafo
        transformer = oms.MapAlignmentTransformer()
        transformer.transformRetentionTimes(feature_map, trafo, True)

    return trafos

# align mzML files based on FeatureMap alignment and store as mzML files (for GNPS!)
def generateAlignedMzML(trafos, mzML_files):
    for file in mzML_files:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(file, exp)
        exp.sortSpectra(True)
        exp.setMetaValue("mzML_path", file)
        if file not in trafos.keys():
            oms.MzMLFile().store(file[:-5] + "_aligned.mzML", exp)
            continue
        transformer = oms.MapAlignmentTransformer()
        trafo_description = trafos[file]
        transformer.transformRetentionTimes(exp, trafo_description, True)
        oms.MzMLFile().store(file[:-5] + "_aligned.mzML", exp)
    mzML_files = [file[:-5] + "_aligned.mzML" for file in mzML_files]
    return mzML_files

def mapMS2(mzML_files, feature_maps):
    feature_maps_mapped = []
    use_centroid_rt = False
    use_centroid_mz = True
    mapper = oms.IDMapper()
    for file in mzML_files:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(file, exp)
        for i, feature_map in enumerate(feature_maps):
            if feature_map.getMetaValue("spectra_data")[
                0
            ].decode() == exp.getMetaValue("mzML_path"):
                peptide_ids = []
                protein_ids = []
                mapper.annotate(
                    feature_map,
                    peptide_ids,
                    protein_ids,
                    use_centroid_rt,
                    use_centroid_mz,
                    exp,
                )
                fm_new = oms.FeatureMap(feature_map)
                fm_new.clear(False)
                # set unique identifiers to protein and peptide identifications
                prot_ids = []
                if len(feature_map.getProteinIdentifications()) > 0:
                    prot_id = feature_map.getProteinIdentifications()[0]
                    prot_id.setIdentifier(f"Identifier_{i}")
                    prot_ids.append(prot_id)
                fm_new.setProteinIdentifications(prot_ids)
                for feature in feature_map:
                    pep_ids = []
                    for pep_id in feature.getPeptideIdentifications():
                        pep_id.setIdentifier(f"Identifier_{i}")
                        pep_ids.append(pep_id)
                    feature.setPeptideIdentifications(pep_ids)
                    fm_new.push_back(feature)
                feature_maps_mapped.append(fm_new)
    feature_maps = feature_maps_mapped
    return feature_maps

def getAdducts(feature_maps):
    feature_maps_adducts = []
    for feature_map in feature_maps:
        mfd = oms.MetaboliteFeatureDeconvolution()
        mdf_par = mfd.getDefaults()
        mdf_par.setValue(
            "potential_adducts",
            [
                b"H:+:0.4",
                b"Na:+:0.2",
                b"NH4:+:0.2",
                b"H-1O-1:+:0.1",
                b"H-3O-2:+:0.1",
            ],
        )
        mfd.setParameters(mdf_par)
        feature_map_adduct = oms.FeatureMap()
        mfd.compute(feature_map, feature_map_adduct, oms.ConsensusMap(), oms.ConsensusMap())
        feature_maps_adducts.append(feature_map_adduct)
    feature_maps = feature_maps_adducts
    
    # for SIRIUS store the feature maps as featureXML files!
    for feature_map in feature_maps:
        oms.FeatureXMLFile().store(
            feature_map.getMetaValue("spectra_data")[0].decode()[:-4]
            + "featureXML",
            feature_map,
        )
    return feature_maps

def generateConsensusMap(feature_maps, consensus_filename="FeatureMatrix.consensusXML"):
    feature_grouper = oms.FeatureGroupingAlgorithmKD()
    
    consensus_map = oms.ConsensusMap()
    file_descriptions = consensus_map.getColumnHeaders()
    
    for i, feature_map in enumerate(feature_maps):
        file_description = file_descriptions.get(i, oms.ColumnHeader())
        file_description.filename = os.path.basename(
            feature_map.getMetaValue("spectra_data")[0].decode()
        )
        file_description.size = feature_map.size()
        file_descriptions[i] = file_description
    
    feature_grouper.group(feature_maps, consensus_map)
    consensus_map.setColumnHeaders(file_descriptions)
    consensus_map.setUniqueIds()
    oms.ConsensusXMLFile().store(consensus_filename, consensus_map)
    return consensus_map

def filterConsensus(input_consensus='FeatureMatrix.consensusXML', out_consensus="filtered.consensusXML"):
    consensusXML_file = input_consensus
    
    consensus_map = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(consensusXML_file, consensus_map)
    filtered_map = oms.ConsensusMap(consensus_map)
    filtered_map.clear(False)
    for feature in consensus_map:
        if feature.getPeptideIdentifications():
            filtered_map.push_back(feature)
    
    consensusXML_file = out_consensus
    oms.ConsensusXMLFile().store(consensusXML_file, filtered_map)

def export2GNPS(mzML_files, consensus_map, consensusXML_file, out_dir='.'):
    oms.GNPSMGFFile().store(
        oms.String(consensusXML_file),
        [file.encode() for file in mzML_files],
        oms.String("MS2data.mgf"),
    )
    oms.GNPSQuantificationFile().store(consensus_map, f"{out_dir}/FeatureQuantificationTable.txt")
    oms.GNPSMetaValueFile().store(consensus_map, f"{out_dir}/MetaValueTable.tsv")

def export2IIMN(consensus_map, out_dir='.'):
    # for IIMN
    oms.IonIdentityMolecularNetworking().annotateConsensusMap(consensus_map)
    oms.IonIdentityMolecularNetworking().writeSupplementaryPairTable(
        consensus_map, f"{out_dir}/SupplementaryPairTable.csv"
    )


