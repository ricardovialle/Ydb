library(shiny);library(shinyjs);library(markdown);library(DT);library(ggplot2);library(scales)
library(rCharts);library(shinythemes);library(shinyWidgets);library(shinydashboard)
library(shinydashboardPlus);library(shinycssloaders);library(billboarder);library(plotly)
library(doBy);library(plyr);library(dplyr);library(shinydashboardPlus);library(DBI);library(heatmaply)
library(pool);library(datasets);library(reshape2);library(visNetwork);library(d3heatmap);library(ggsci)
library(RMySQL)

#################################################################################
#### Definitions ################################################################
#################################################################################

#### Atom interaction ####
atom_interaction_type = data.frame(
  code = c("F","E","A","P","V","C","H","I","J","K","L","M","N","O","Q","R","S","T"),
  type = c("Hydrophobic","Electrostatic","Aromatic","Pi-cation",
           "Van der Waals","Van der Waals clash","Hydrogen bond-MM",
           "Hydrogen bond-SM","Hydrogen bond-MS","Hydrogen bond-SS",
           "Hydrogen bond-Water mediated-MM","Hydrogen bond-Water mediated-SM",
           "Hydrogen bond-Water mediated-MS","Hydrogen bond-Water mediated-SS",
           "Hydrogen bond-2 Waters mediated-MM","Hydrogen bond-2 Waters mediated-SM",
           "Hydrogen bond-2 Waters mediated-MS","Hydrogen bond-2 Waters mediated-SS")
)
#### ####

#### AA classes #####
aa_classes = data.frame(
  class = c("Acid","Acid",
            "Basic","Basic","Basic",
            "Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic","Hydrophobic",
            "Neutral","Neutral",
            "Polar","Polar","Polar","Polar","Polar",
            "Gap"
  ), 
  class_color = c("red","red",
                  "blue","blue","blue",
                  "black","black","black","black","black","black","black","black",
                  "purple","purple",
                  "green","green","green","green","green",
                  "gray"
  ),
  aa = c("D","E",
         "R","H","K",
         "F","P","V","L","I","A","W","M",
         "N","Q",
         "C","S","G","Y","T",
         "-")
)
#### ####

#### Complex table - Columns dictionary ####
complex_tbl_dic = c(
  "Id","Complex Id","Complex",
  "Pdb","PDB","Complex",
  "Model","Model","Complex",
  "DeltaG","Delta G","Complex",
  "DeltaSAS","Delta SAS","Complex",
  "Affinity","Affinity","Complex",
  "AffinityMethod","Affinity Method","Complex",
  "AffinityTemperature","Affinity Temperature","Complex",
  "AffinityPmid","Affinity PMID","Complex",
  "PDBDescription","PDB Description","Complex",
  "PDBCompound","PDB Compound","Complex",
  "PDBAuthors","PDB Authors","Complex",
  "PDBReference","PDB Reference","Complex",
  "PDBDepositionDate","PDB Deposition Date","Complex",
  "PDBExperiment","PDB Experiment","Complex",
  "PDBResolution","PDB Resolution","Complex",
  "PDBRFree","R Free","Complex",
  "PDBRFactor","R Factor","Complex",
  "HeavyChainId","Chain Id (Heavy)","Chain",
  "HeavyChainName","Chain Name (Heavy)","Chain",
  "HeavyChainType","Chain Type (Heavy)","Chain",
  "HeavyChainSubclass","Chain Subclass (Heavy)","Chain",
  "HeavySpecies","Chain Species (Heavy)","Chain",
  "HeavyLightType","Light Type (Heavy)","Chain",
  "HeavyResSequence","Sequence (Heavy)","Chain",
  "HeavyAtomSequence","Atom Sequence (Heavy)","Chain",
  "HeavyGravy","Gravy (Heavy)","Chain",
  "HeavyIsoelectricPoint","Isoelectric Point (Heavy)","Chain",
  "HeavyNumberOfAas","Number of Amino Acids (Heavy)","Chain",
  "HeavyTotalNumberofAas","Total Number of AAs (Heavy)","Chain",
  "HeavyAliphaticIndex","Aliphatic Index (Heavy)","Chain",
  "HeavyMolecularWeight","Molecular Weight (Heavy)","Chain",
  "HeavyNegativeCharged","Negative Charged (Heavy)","Chain",
  "HeavyPositiveCharged","Positive Charged (Heavy)","Chain",
  "HeavyPolarUncharged","Polar Uncharged (Heavy)","Chain",
  "HeavySmall","Small (Heavy)","Chain",
  "HeavyHydrophobic","Hydrophobic (Heavy)","Chain",
  "HeavyAlcohol","Alcohol (Heavy)","Chain",
  "HeavyAromatic","Aromatic (Heavy)","Chain",
  "HeavyAlphaHelixStructure","Alpha-Helix (Heavy)","SS",
  "HeavyIsolatedBetaBridgeStructure","Isolated Beta Bridge (Heavy)","SS",
  "HeavyStrandStructure","Strand (Heavy)","SS",
  "HeavyHelix3_10Structure","Helix 3_10 (Heavy)","SS",
  "HeavyPiHelixStructure","Pi Helix (Heavy)","SS",
  "HeavyTurnStructure","Turn (Heavy)","SS",
  "HeavyBendStructure","Bend (Heavy)","SS",
  "HeavyNoneStructure","None secondary structure (Heavy)","SS",
  "HeavyIMGTComment","IMGT Comment (Heavy)","Chain",
  "HeavyPdbComment","PDB Comment (Heavy)","Chain",
  "LightChainId","Chain Id (Light)","Chain",
  "LightChainName","Chain Name (Light)","Chain",
  "LightChainType","Chain Type (Light)","Chain",
  "LightChainSubclass","Chain Subclass (Light)","Chain",
  "LightSpecies","Chain Species (Light)","Chain",
  "LightLightType","Light Type (Light)","Chain",
  "LightResSequence","Sequence (Light)","Chain",
  "LightAtomSequence","Atom Sequence (Light)","Chain",
  "LightGravy","Gravy (Light)","Chain",
  "LightIsoelectricPoint","Isoelectric Point (Light)","Chain",
  "LightNumberOfAas","Number of Amino Acids (Light)","Chain",
  "LightTotalNumberofAas","Total Number of AAs (Light)","Chain",
  "LightAliphaticIndex","Aliphatic Index (Light)","Chain",
  "LightMolecularWeight","Molecular Weight (Light)","Chain",
  "LightNegativeCharged","Negative Charged (Light)","Chain",
  "LightPositiveCharged","Positive Charged (Light)","Chain",
  "LightPolarUncharged","Polar Uncharged (Light)","Chain",
  "LightSmall","Small (Light)","Chain",
  "LightHydrophobic","Hydrophobic (Light)","Chain",
  "LightAlcohol","Alcohol (Light)","Chain",
  "LightAromatic","Aromatic (Light)","Chain",
  "LightAlphaHelixStructure","Alpha-Helix (Light)","SS",
  "LightIsolatedBetaBridgeStructure","Isolated Beta Bridge (Light)","SS",
  "LightStrandStructure","Strand (Light)","SS",
  "LightHelix3_10Structure","Helix 3_10 (Light)","SS",
  "LightPiHelixStructure","Pi Helix (Light)","SS",
  "LightTurnStructure","Turn (Light)","SS",
  "LightBendStructure","Bend (Light)","SS",
  "LightNoneStructure","None secondary structure (Light)","SS",
  "LightIMGTComment","IMGT Comment (Light)","Chain",
  "LightPdbComment","PDB Comment (Light)","Chain",
  "AntigenName","Antigen Name","Chain",
  "AntigenChainId","Chain Id (Antigen)","Chain",
  "AntigenChainName","Chain Name (Antigen)","Chain",
  "AntigenChainType","Chain Type (Antigen)","Chain",
  "AntigenChainSubclass","Chain Subclass (Antigen)","Chain",
  "AntigenSpecies","Chain Species (Antigen)","Chain",
  "AntigenLightType","Light Type (Antigen)","Chain",
  "AntigenResSequence","Sequence (Antigen)","Chain",
  "AntigenAtomSequence","Atom Sequence (Antigen)","Chain",
  "AntigenGravy","Gravy (Antigen)","Chain",
  "AntigenIsoelectricPoint","Isoelectric Point (Antigen)","Chain",
  "AntigenNumberOfAas","Number of Amino Acids (Antigen)","Chain",
  "AntigenTotalNumberofAas","Total Number of AAs (Antigen)","Chain",
  "AntigenAliphaticIndex","Aliphatic Index (Antigen)","Chain",
  "AntigenMolecularWeight","Molecular Weight (Antigen)","Chain",
  "AntigenNegativeCharged","Negative Charged (Antigen)","Chain",
  "AntigenPositiveCharged","Positive Charged (Antigen)","Chain",
  "AntigenPolarUncharged","Polar Uncharged (Antigen)","Chain",
  "AntigenSmall","Small (Antigen)","Chain",
  "AntigenHydrophobic","Hydrophobic (Antigen)","Chain",
  "AntigenAlcohol","Alcohol (Antigen)","Chain",
  "AntigenAromatic","Aromatic (Antigen)","Chain",
  "AntigenAlphaHelixStructure","Alpha-Helix (Antigen)","SS",
  "AntigenIsolatedBetaBridgeStructure","Isolated Beta Bridge (Antigen)","SS",
  "AntigenStrandStructure","Strand (Antigen)","SS",
  "AntigenHelix3_10Structure","Helix 3_10 (Antigen)","SS",
  "AntigenPiHelixStructure","Pi Helix (Antigen)","SS",
  "AntigenTurnStructure","Turn (Antigen)","SS",
  "AntigenBendStructure","Bend (Antigen)","SS",
  "AntigenNoneStructure","None secondary structure (Antigen)","SS",
  "AntigenIMGTComment","IMGT Comment (Antigen)","Chain",
  "AntigenPdbComment","PDB Comment (Antigen)","Chain",
  "Interaction","Interaction","Complex"
)
complex_tbl_dic = data.frame(colname = complex_tbl_dic[seq(1,length(complex_tbl_dic),by = 3)],
                             desc = complex_tbl_dic[seq(2,length(complex_tbl_dic),by = 3)],
                             class = complex_tbl_dic[seq(3,length(complex_tbl_dic),by = 3)])
#### ####

#### Chain table - Columns dictionary ####
chain_tbl_dic = c(
  "Id","Chain Id","Chain",
  "Pdb","PDB","Chain",
  "ChainName","Chain Name","Chain",
  "ChainType","Chain Type","Chain",
  "SpeciesId","Species Id","Chain",
  "Species","Species","Chain",
  "Gravy","Gravy","Chain",
  "IsoelectricPoint","Isoelectric Point","Chain",
  "AliphaticIndex","Aliphatic Index","Chain",
  "MolecularWeight","Molecular Weight","Chain",
  "TotalNumberOfAas","Number of AAs (Total)","Chain",
  "NumberOfAas","Number of AAs","Chain",
  "NegativeCharged","Negative Charged","Chain",
  "PositiveCharged","Positive Charged","Chain",
  "PolarUncharged","Polar Uncharged","Chain",
  "Small","Small","Chain",
  "Hydrophobic","Hydrophobic","Chain",
  "Alcohol","Alcohol","Chain",
  "Aromatic","Aromatic","Chain",
  "AlphaHelixStructure","Alpha-Helix","SS",        
  "IsolatedBetaBridgeStructure","Isolated Beta Bridge","SS",
  "StrandStructure","Strand","SS",
  "Helix3_10Structure","Helix 3-10","SS",
  "PiHelixStructure","Pi Helix","SS",
  "TurnStructure","Turn","SS",
  "BendStructure","Bend","SS",
  "NoneStructure","None secondary structure","SS"        
)
chain_tbl_dic = data.frame(colname = chain_tbl_dic[seq(1,length(chain_tbl_dic),by = 3)],
                           desc = chain_tbl_dic[seq(2,length(chain_tbl_dic),by = 3)],
                           class = chain_tbl_dic[seq(3,length(chain_tbl_dic),by = 3)])
#### ####

#### Interaction table - Columns dictionary ####
interaction_tbl_dic = c(
  "AtomInteraction_AntibodyAtom","Antibody Atom","Antibody",
  "AtomInteraction_AntibodyResidueId","Antibody Residue","Antibody",
  "AntibodyResidue_Id","Antibody Res Id","Antibody",
  "AntibodyResidue_ChainId","Antibody Chain Id","Antibody",
  "AntibodyResidue_PdbIndex","Residue PDB Index (Antibody)","Antibody",
  "AntibodyResidue_PdbIndexCode","Residue PDB Index code (Antibody)","Antibody",
  "AntibodyResidue_AminoAcid","Amino acid (Antibody)","Antibody",
  "AntibodyResidue_Flexibility","Residue Flexibility (Antibody)","Antibody",
  "AntibodyResidue_FlexibilityConfidence","Flexibility Confidence (Antibody)","Antibody",
  "AntibodyResidue_StructuralProperty","Residue Structural Property (Antibody)","Antibody",
  "AntibodyResidue_SASComplexed","Residue SAS Complexed (Antibody)","Antibody",
  "AntibodyResidue_SASUncomplexed","Residue SAS Uncomplexed (Antibody)","Antibody",
  "AntibodyResidue_IMGTIndex","Residue IMGT Index (Antibody)","Antibody",
  "AtomInteraction_AntigenAtom","Antigen Atom","Antigen",
  "AtomInteraction_AntigenResidueId","Antigen Residue","Antigen",
  "AntigenResidue_Id","Antigen Res Id","Antigen",
  "AntigenResidue_ChainId","Antigen Chain Id","Antigen",
  "AntigenResidue_PdbIndex","Residue PDB Index (Antigen)","Antigen",
  "AntigenResidue_PdbIndexCode","Residue PDB Index code (Antigen)","Antigen",
  "AntigenResidue_AminoAcid","Amino acid (Antigen)","Antigen",
  "AntigenResidue_Flexibility","Residue Flexibility (Antigen)","Antigen",
  "AntigenResidue_FlexibilityConfidence","Flexibility Confidence (Antigen)","Antigen",
  "AntigenResidue_StructuralProperty","Residue Structural Property (Antigen)","Antigen",
  "AntigenResidue_SASComplexed","Residue SAS Complexed (Antigen)","Antigen",
  "AntigenResidue_SASUncomplexed","Residue SAS Uncomplexed (Antigen)","Antigen",
  "AntigenResidue_IMGTIndex","Residue IMGT Index (Antigen)","Antigen",
  "AtomInteraction_Type","Interaction Type","Interaction",
  "AtomInteraction_Distance","Interaction Distance","Interaction"
)
interaction_tbl_dic = data.frame(colname = interaction_tbl_dic[seq(1,length(interaction_tbl_dic),by = 3)],
                                 desc = interaction_tbl_dic[seq(2,length(interaction_tbl_dic),by = 3)],
                                 class = interaction_tbl_dic[seq(3,length(interaction_tbl_dic),by = 3)])
#### ####

#### CATH table ####
#cathdb = read.table("data/YdbCATHDomainResidue.txt", sep = '\t')
#colnames(cathdb) = c("PdbId","ChainId","Domain","ChathId","Class","Architecture","Topology","Homologous","Desc_Class","Desc_Arch","Desc_Topo","Desc_Homol","BeginPos","EndPos")
#### ####

#################################################################################
#### MySQL Connection ###########################################################
#################################################################################

#### Connection pool ####
pool <- dbPool(
  drv = RMySQL::MySQL(),
  dbname = "Ydb",
  host = "127.0.0.1",
  username = "ydbInterface",
  password = "senhaYdb")

killDbConnections <- function () {
  all_cons <- dbListConnections(RMySQL::MySQL())
  print(all_cons)
  for(con in all_cons)
    +  dbDisconnect(con)
  print(paste(length(all_cons), " connections killed."))
}

onStop(function() {
  poolClose(pool)
})
#### #### 

#################################################################################
#### MySQL Queries ##############################################################
#################################################################################

#### Function: Complex table (all) ####
complex_table_show_all <- function(){
  sql <- "SELECT C.Id, C.Pdb, C.Model, C.DeltaG, C.DeltaSAS, C.Affinity, C.AffinityMethod, C.AffinityTemperature, C.AffinityPmid,
  
  PE.Description AS PDBDescription, PE.Compound AS PDBCompound,
  PE.Authors AS PDBAuthors, PE.Reference AS PDBReference, PE.DepositionDate AS PDBDepositionDate,
  PE.Experiment AS PDBExperiment, PE.Resolution AS PDBResolution,
  PE.RFree AS PDBRFree, PE.RFactor AS PDBRFactor,
  
  HC.Id AS HeavyChainId,
  HC.ChainName AS HeavyChainName, HC.ChainType AS HeavyChainType, 
  HC.ChainSubclass AS HeavyChainSubclass, SHC.Name AS HeavySpecies, HC.LightType AS HeavyLightType, 
  HC.ResSequence AS HeavyResSequence, HC.AtomSequence AS HeavyAtomSequence, 
  HC.Gravy AS HeavyGravy, HC.IsoelectricPoint AS HeavyIsoelectricPoint, HC.NumberOfAas AS HeavyNumberOfAas,
  HC.TotalNumberofAas AS HeavyTotalNumberofAas, HC.AliphaticIndex AS HeavyAliphaticIndex,
  HC.MolecularWeight AS HeavyMolecularWeight, HC.NegativeCharged AS HeavyNegativeCharged,
  HC.PositiveCharged AS HeavyPositiveCharged, HC.PolarUncharged AS HeavyPolarUncharged, HC.Small AS HeavySmall,
  HC.Hydrophobic AS HeavyHydrophobic, HC.Alcohol AS HeavyAlcohol, HC.Aromatic AS HeavyAromatic,
  HC.AlphaHelixStructure AS HeavyAlphaHelixStructure, HC.IsolatedBetaBridgeStructure AS HeavyIsolatedBetaBridgeStructure,
  HC.StrandStructure AS HeavyStrandStructure, HC.Helix3_10Structure AS HeavyHelix3_10Structure, HC.PiHelixStructure AS HeavyPiHelixStructure,
  HC.TurnStructure AS HeavyTurnStructure, HC.BendStructure AS HeavyBendStructure,              
  HC.NoneStructure AS HeavyNoneStructure, HC.IMGTComment AS HeavyIMGTComment, HC.PdbComment AS HeavyPdbComment, 
  
  LC.Id AS LightChainId,
  LC.ChainName AS LightChainName, LC.ChainType AS LightChainType, 
  LC.ChainSubclass AS LightChainSubclass, SLC.Name AS LightSpecies, LC.LightType AS LightLightType, 
  LC.ResSequence AS LightResSequence, LC.AtomSequence AS LightAtomSequence, 
  LC.Gravy AS LightGravy, LC.IsoelectricPoint AS LightIsoelectricPoint, LC.NumberOfAas AS LightNumberOfAas,
  LC.TotalNumberofAas AS LightTotalNumberofAas, LC.AliphaticIndex AS LightAliphaticIndex,
  LC.MolecularWeight AS LightMolecularWeight, LC.NegativeCharged AS LightNegativeCharged,
  LC.PositiveCharged AS LightPositiveCharged, LC.PolarUncharged AS LightPolarUncharged, LC.Small AS LightSmall,
  LC.Hydrophobic AS LightHydrophobic, LC.Alcohol AS LightAlcohol, LC.Aromatic AS LightAromatic,
  LC.AlphaHelixStructure AS LightAlphaHelixStructure, LC.IsolatedBetaBridgeStructure AS LightIsolatedBetaBridgeStructure,
  LC.StrandStructure AS LightStrandStructure, LC.Helix3_10Structure AS LightHelix3_10Structure, LC.PiHelixStructure AS LightPiHelixStructure,
  LC.TurnStructure AS LightTurnStructure, LC.BendStructure AS LightBendStructure,              
  LC.NoneStructure AS LightNoneStructure, LC.IMGTComment AS LightIMGTComment, LC.PdbComment AS LightPdbComment,
  
  
  A.AntigenName AS AntigenName,
  AC.Id AS AntigenChainId,
  AC.ChainName AS AntigenChainName, AC.ChainType AS AntigenChainType, 
  AC.ChainSubclass AS AntigenChainSubclass, SAC.Name AS AntigenSpecies, AC.LightType AS AntigenLightType, 
  AC.ResSequence AS AntigenResSequence, AC.AtomSequence AS AntigenAtomSequence, 
  AC.Gravy AS AntigenGravy, AC.IsoelectricPoint AS AntigenIsoelectricPoint, AC.NumberOfAas AS AntigenNumberOfAas,
  AC.TotalNumberofAas AS AntigenTotalNumberofAas, AC.AliphaticIndex AS AntigenAliphaticIndex,
  AC.MolecularWeight AS AntigenMolecularWeight, AC.NegativeCharged AS AntigenNegativeCharged,
  AC.PositiveCharged AS AntigenPositiveCharged, AC.PolarUncharged AS AntigenPolarUncharged, AC.Small AS AntigenSmall,
  AC.Hydrophobic AS AntigenHydrophobic, AC.Alcohol AS AntigenAlcohol, AC.Aromatic AS AntigenAromatic,
  AC.AlphaHelixStructure AS AntigenAlphaHelixStructure, AC.IsolatedBetaBridgeStructure AS AntigenIsolatedBetaBridgeStructure,
  AC.StrandStructure AS AntigenStrandStructure, AC.Helix3_10Structure AS AntigenHelix3_10Structure, AC.PiHelixStructure AS AntigenPiHelixStructure,
  AC.TurnStructure AS AntigenTurnStructure, AC.BendStructure AS AntigenBendStructure,              
  AC.NoneStructure AS AntigenNoneStructure, AC.IMGTComment AS AntigenIMGTComment, AC.PdbComment AS AntigenPdbComment,
  
  'Show details' AS Interaction
  
  FROM Complex AS C 
  LEFT JOIN Chain AS HC ON C.HeavyChainId = HC.Id 
  LEFT JOIN Species AS SHC ON SHC.Id = HC.SpeciesId
  LEFT JOIN Chain AS LC ON C.LightChainId = LC.Id
  LEFT JOIN Species AS SLC ON SLC.Id = LC.SpeciesId
  LEFT JOIN AntigenChain AS A ON C.Id = A.ComplexId
  LEFT JOIN Chain AS AC ON A.ChainId = AC.Id 
  LEFT JOIN Species AS SAC ON SAC.Id = AC.SpeciesId
  LEFT JOIN PdbEntry AS PE ON C.Pdb = PE.Pdb 
  WHERE ?id ORDER BY C.Id ASC;"
  
  query <- sqlInterpolate(pool, sql, id = 1)
  complex_table = dbGetQuery(pool, query)
  return(complex_table)
}
#### #### 

#### Function: Complex table (selected) ####
complex_table_show_selected <- function(params){
  sql <- paste0("SELECT C.Id, C.Pdb, C.Model, C.DeltaG, C.DeltaSAS, C.Affinity, C.AffinityMethod, C.AffinityTemperature, C.AffinityPmid,
  
  PE.Description AS PDBDescription, PE.Compound AS PDBCompound,
  PE.Authors AS PDBAuthors, PE.Reference AS PDBReference, PE.DepositionDate AS PDBDepositionDate,
  PE.Experiment AS PDBExperiment, PE.Resolution AS PDBResolution,
  PE.RFree AS PDBRFree, PE.RFactor AS PDBRFactor,
  
  HC.Id AS HeavyChainId,
  HC.ChainName AS HeavyChainName, HC.ChainType AS HeavyChainType, 
  HC.ChainSubclass AS HeavyChainSubclass, SHC.Name AS HeavySpecies, HC.LightType AS HeavyLightType, 
  HC.ResSequence AS HeavyResSequence, HC.AtomSequence AS HeavyAtomSequence, 
  HC.Gravy AS HeavyGravy, HC.IsoelectricPoint AS HeavyIsoelectricPoint, HC.NumberOfAas AS HeavyNumberOfAas,
  HC.TotalNumberofAas AS HeavyTotalNumberofAas, HC.AliphaticIndex AS HeavyAliphaticIndex,
  HC.MolecularWeight AS HeavyMolecularWeight, HC.NegativeCharged AS HeavyNegativeCharged,
  HC.PositiveCharged AS HeavyPositiveCharged, HC.PolarUncharged AS HeavyPolarUncharged, HC.Small AS HeavySmall,
  HC.Hydrophobic AS HeavyHydrophobic, HC.Alcohol AS HeavyAlcohol, HC.Aromatic AS HeavyAromatic,
  HC.AlphaHelixStructure AS HeavyAlphaHelixStructure, HC.IsolatedBetaBridgeStructure AS HeavyIsolatedBetaBridgeStructure,
  HC.StrandStructure AS HeavyStrandStructure, HC.Helix3_10Structure AS HeavyHelix3_10Structure, HC.PiHelixStructure AS HeavyPiHelixStructure,
  HC.TurnStructure AS HeavyTurnStructure, HC.BendStructure AS HeavyBendStructure,              
  HC.NoneStructure AS HeavyNoneStructure, HC.IMGTComment AS HeavyIMGTComment, HC.PdbComment AS HeavyPdbComment, 
  
  LC.Id AS LightChainId,
  LC.ChainName AS LightChainName, LC.ChainType AS LightChainType, 
  LC.ChainSubclass AS LightChainSubclass, SLC.Name AS LightSpecies, LC.LightType AS LightLightType, 
  LC.ResSequence AS LightResSequence, LC.AtomSequence AS LightAtomSequence, 
  LC.Gravy AS LightGravy, LC.IsoelectricPoint AS LightIsoelectricPoint, LC.NumberOfAas AS LightNumberOfAas,
  LC.TotalNumberofAas AS LightTotalNumberofAas, LC.AliphaticIndex AS LightAliphaticIndex,
  LC.MolecularWeight AS LightMolecularWeight, LC.NegativeCharged AS LightNegativeCharged,
  LC.PositiveCharged AS LightPositiveCharged, LC.PolarUncharged AS LightPolarUncharged, LC.Small AS LightSmall,
  LC.Hydrophobic AS LightHydrophobic, LC.Alcohol AS LightAlcohol, LC.Aromatic AS LightAromatic,
  LC.AlphaHelixStructure AS LightAlphaHelixStructure, LC.IsolatedBetaBridgeStructure AS LightIsolatedBetaBridgeStructure,
  LC.StrandStructure AS LightStrandStructure, LC.Helix3_10Structure AS LightHelix3_10Structure, LC.PiHelixStructure AS LightPiHelixStructure,
  LC.TurnStructure AS LightTurnStructure, LC.BendStructure AS LightBendStructure,              
  LC.NoneStructure AS LightNoneStructure, LC.IMGTComment AS LightIMGTComment, LC.PdbComment AS LightPdbComment,
  
  
  A.AntigenName AS AntigenName,
  AC.Id AS AntigenChainId,
  AC.ChainName AS AntigenChainName, AC.ChainType AS AntigenChainType, 
  AC.ChainSubclass AS AntigenChainSubclass, SAC.Name AS AntigenSpecies, AC.LightType AS AntigenLightType, 
  AC.ResSequence AS AntigenResSequence, AC.AtomSequence AS AntigenAtomSequence, 
  AC.Gravy AS AntigenGravy, AC.IsoelectricPoint AS AntigenIsoelectricPoint, AC.NumberOfAas AS AntigenNumberOfAas,
  AC.TotalNumberofAas AS AntigenTotalNumberofAas, AC.AliphaticIndex AS AntigenAliphaticIndex,
  AC.MolecularWeight AS AntigenMolecularWeight, AC.NegativeCharged AS AntigenNegativeCharged,
  AC.PositiveCharged AS AntigenPositiveCharged, AC.PolarUncharged AS AntigenPolarUncharged, AC.Small AS AntigenSmall,
  AC.Hydrophobic AS AntigenHydrophobic, AC.Alcohol AS AntigenAlcohol, AC.Aromatic AS AntigenAromatic,
  AC.AlphaHelixStructure AS AntigenAlphaHelixStructure, AC.IsolatedBetaBridgeStructure AS AntigenIsolatedBetaBridgeStructure,
  AC.StrandStructure AS AntigenStrandStructure, AC.Helix3_10Structure AS AntigenHelix3_10Structure, AC.PiHelixStructure AS AntigenPiHelixStructure,
  AC.TurnStructure AS AntigenTurnStructure, AC.BendStructure AS AntigenBendStructure,              
  AC.NoneStructure AS AntigenNoneStructure, AC.IMGTComment AS AntigenIMGTComment, AC.PdbComment AS AntigenPdbComment,
  
  'Show details' AS Interaction
  
  FROM Complex AS C 
  LEFT JOIN Chain AS HC ON C.HeavyChainId = HC.Id 
  LEFT JOIN Species AS SHC ON SHC.Id = HC.SpeciesId
  LEFT JOIN Chain AS LC ON C.LightChainId = LC.Id
  LEFT JOIN Species AS SLC ON SLC.Id = LC.SpeciesId
  LEFT JOIN AntigenChain AS A ON C.Id = A.ComplexId
  LEFT JOIN Chain AS AC ON A.ChainId = AC.Id 
  LEFT JOIN Species AS SAC ON SAC.Id = AC.SpeciesId
  LEFT JOIN PdbEntry AS PE ON C.Pdb = PE.Pdb 
  WHERE ?id ORDER BY C.Id ASC;")
  
  query <- sqlInterpolate(pool, sql, id = 1)
  complex_table = dbGetQuery(pool, query)
  return(complex_table)
}
#### #### 

#### Function: Select from table #### 
select_table <- function(engineered_antibodies,
                         right_amino_acid_in_conserved_positions,
                         only_one_VHVL_pair,
                         identity_filter,
                         identity_cutoff){
  return(paste0(engineered_antibodies,
                right_amino_acid_in_conserved_positions,
                only_one_VHVL_pair,
                identity_filter,
                identity_cutoff,
                collapse = "|"))
}
#### #### 

#################################################################################
#### Define server logic ########################################################
#################################################################################
shinyServer(function(input, output, session) {

  observeEvent(input$actionButton.all_pdb_checkbox, {
    # Check if at least one opions is selected. Show message if not.
    if (length(input$all_pdb_checkbox)==0){
      session$sendCustomMessage(type = 'testmessage',
                                message = 'Select an option.')
      return()
    }
    # Get checkbox values
    engineered_antibodies = input$all_pdb_checkbox[1]
    right_amino_acid_in_conserved_positions = input$all_pdb_checkbox[2]
    only_one_VHVL_pair = input$all_pdb_checkbox[3]
    identity_filter = input$all_pdb_checkbox[4]
    identity_cutoff = ifelse(!is.na(identity_filter),input$identity_cutoff,0)
    # Change to table page
    updateTabItems(session, "menu", "complex")
    # Return options to interface
    output$all_pdb_checkbox_sel <- renderPrint({ 
      return(select_table(engineered_antibodies,
                          right_amino_acid_in_conserved_positions,
                          only_one_VHVL_pair,
                          identity_filter,
                          identity_cutoff))
    })
  })

  #### Render the table containing shiny inputs ####
  output$complex_table <- DT::renderDataTable({
    # Get current data 
    data <- mtcars
    # The table
    datatable(data, 
              style = 'bootstrap',
              rownames= FALSE,
              selection = list(mode = "multiple", target= 'row'),
              options = list(
                pageLength = 10,
                scrollX = TRUE,
                deferRender = TRUE
              )
    )
  }, server = F)
  #### ####
  
})
