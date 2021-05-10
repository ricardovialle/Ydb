fname_list = "cath-domain-list.txt"
tbl_list = read.table(fname_list, comment.char = '#', sep = "", quote = "")
colnames(tbl_list) = c("cathId","class","architecture","topology","homologous","S35","S60","S95","S100","S100_count","domain_len","resolution")

fname_names = "cath-names.txt"
cathNames = data.frame()
tbl_names = readLines(fname_names)
k = 0
for (i in 1:length(tbl_names)) {
   if (!substr(tbl_names[i],1,1)=='#') {
      line = tbl_names[i]
      desc_split = strsplit(line,':')[[1]]
      desc = desc_split[2]
      node_split = strsplit(desc_split[1],'\\s+')[[1]]
      nodeNum = node_split[1]
      nodeNum_split = strsplit(nodeNum,'\\.')[[1]]
      k = k + 1
      cathNames[k,"nodeNum"] = nodeNum
      cathNames[k,"class"] = nodeNum_split[1]
      cathNames[k,"architecture"] = nodeNum_split[2]
      cathNames[k,"topology"] = nodeNum_split[3]
      cathNames[k,"homologous"] = nodeNum_split[4]
      cathNames[k,"Desc"] = desc
   }
}

fname = "cath-domain-boundaries-seqreschopping.txt"
tbl = read.table(fname, comment.char = '#', sep = '\t', quote = "")
tbl_split = transform(tbl, PDB = substr(V1,1,4), Chain = substr(V1,5,5), Domain = substr(V1,6,7))
for (i in 1:nrow(tbl_split)) {
   t = as.character(tbl_split[i,"V2"])
   if (grepl(",", t)){
      dom_split = strsplit(t,',')[[1]]
   }else{
      dom_split = t
   }
   print(i)
   for (j in 1:length(dom_split)){
      dom_pos = strsplit(dom_split[j],'-')[[1]]
      pbegin = as.numeric(dom_pos[1])
      pend = as.numeric(dom_pos[2])
      CathId = as.character(tbl_split$V1[1])
      cathcodes = tbl_list[which(as.character(tbl_list$cathId)==CathId),]
      DomainResidue = data.frame(PdbId = tbl_split$PDB[i],
         ChainId = as.character(tbl_split$Chain[i]),
         Domain = as.character(tbl_split$Domain[i]),
         CathId = CathId,
         Class = cathcodes$class,
         Architecture = cathcodes$architecture,
         Topology = cathcodes$topology,
         Homologous = cathcodes$homologous, 
         Desc_Class = cathNames[which(cathNames$class==cathcodes$class & is.na(cathNames$architecture) & is.na(cathNames$topology) & is.na(cathNames$homologous)),"Desc"],
         Desc_Arch = cathNames[which(cathNames$class==cathcodes$class & cathNames$architecture==cathcodes$architecture & is.na(cathNames$topology) & is.na(cathNames$homologous)),"Desc"],
         Desc_Topology = cathNames[which(cathNames$class==cathcodes$class & cathNames$architecture==cathcodes$architecture & cathNames$topology==cathcodes$topology & is.na(cathNames$homologous)),"Desc"],
         Desc_Homologous = cathNames[which(cathNames$class==cathcodes$class & cathNames$architecture==cathcodes$architecture & cathNames$topology==cathcodes$topology & cathNames$homologous==cathcodes$homologous),"Desc"],
         BeginPos = pbegin,
         EndPos = pend)
      write.table(DomainResidue, 
         file = "YdbCATHDomainResidue.txt", sep = '\t',
         quote = F, row.names = F, col.names = F,
         append = T)
   }
}





