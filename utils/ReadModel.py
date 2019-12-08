import math

def ReadFONLL(FileNameFONLL) :
  Col = [[] for iCol in range(16)]
  with open(FileNameFONLL) as FileFONLL:
    for iLine in FileFONLL:
      if "#" not in iLine:
        iLine.strip().split(',')
        iLineSep = iLine.split()
        if len(iLineSep) > 15:
          for iCol in range(16) :
            Col[iCol].append(float(iLineSep[iCol]))
  FONLLdic = {'PtMin' : Col[0],'Cent' : Col[1],'Min' : Col[2],'Max' : Col[3],'Min_sc' : Col[4],'Max_sc' : Col[5],'Min_Mass' : Col[6],'Max_Mass' : Col[7],'Min_PDF' : Col[8],'Max_PDF' : Col[9],'Fr_dot5dot5' : Col[10],'Fr_22' : Col[11],'Fr_21' : Col[12],'Fr_12' : Col[13],'Fr_1dot5' : Col[14],'Fr_dot51' : Col[15]}
  return FONLLdic

def ReadGMVFNS(FileNameGMVFNS, isSACOT=False) :
  Col = [[] for iCol in range(4)]
  with open(FileNameGMVFNS) as FileGMVFNS:
    for iLine in FileGMVFNS:
      if "#" not in iLine:
        iLine.strip().split(',')
        iLineSep = iLine.split()
        if len(iLineSep) > 3:
          if isSACOT:
            Col[0].append(float(iLineSep[0]))
            Col[1].append(float(iLineSep[1]))
            Col[2].append(float(iLineSep[1])-math.sqrt(float(iLineSep[2])**2+float(iLineSep[3])**2))
            Col[3].append(float(iLineSep[1])+math.sqrt(float(iLineSep[2])**2+float(iLineSep[4])**2))
          else:
            for iCol in range(4) :
              Col[iCol].append(float(iLineSep[iCol]))
  if isSACOT:
    GMVFNSdic = {'PtCent' : Col[0],'Cent' : Col[1],'Min' : Col[2],'Max' : Col[3]}
  else:
    GMVFNSdic = {'PtMin' : Col[0],'Cent' : Col[1],'Max' : Col[2],'Min' : Col[3]}
  return GMVFNSdic

def ReadKtFact(FileNameKtFact) :
  Col = [[] for iCol in range(5)]
  with open(FileNameKtFact) as FileKtFact:
    for iLine in FileKtFact:
      if "#" not in iLine:
        iLine.strip().split(',')
        iLineSep = iLine.split()
        if len(iLineSep) > 4:
          for iCol in range(5) :
            Col[iCol].append(float(iLineSep[iCol]))
  KtFactdic = {'PtMin' : Col[0], 'PtMax' : Col[1], 'Cent' : Col[2],'Min' : Col[3],'Max' : Col[4]}
  return KtFactdic

def ReadTAMU(FileNameTAMU) :
  Col = [[] for iCol in range(3)]
  with open(FileNameTAMU) as FileTAMU:
    for iLine in FileTAMU:
      if "#" not in iLine:
        iLine.strip().split(',')
        iLineSep = iLine.split()
        if len(iLineSep) > 2:
          for iCol in range(3) :
            Col[iCol].append(float(iLineSep[iCol]))
  TAMUdic = {'PtCent' : Col[0],'Max' : Col[1],'Min' : Col[2]}
  return TAMUdic

def ReadPHSD(FileNamePHSD) :
  Col = [[] for iCol in range(2)]
  with open(FileNamePHSD) as FilePHSD:
    for iLine in FilePHSD:
      if "#" not in iLine:
        iLine.strip().split(',')
        iLineSep = iLine.split()
        if len(iLineSep) > 1:          
          for iCol in range(2) :
            Col[iCol].append(float(iLineSep[iCol]))
  PHSDdic = {'PtCent' : Col[0],'Cent' : Col[1]}
  return PHSDdic

def ReadGossiaux(FileNameGossiaux) :
  Col = [[] for iCol in range(4)]
  with open(FileNameGossiaux) as FileGossiaux:
    for iLine in FileGossiaux:
      if "#" not in iLine:
        iLine.strip().split(',')
        iLineSep = iLine.split()
        if len(iLineSep) > 3:
          for iCol in range(4) :
            Col[iCol].append(float(iLineSep[iCol]))
  Gossiauxdic = {'PtCent' : Col[0],'Col' : Col[1],'ColRad' : Col[2],'ColRadGluDamp' : Col[3]}
  return Gossiauxdic

def ReadCatania(FileNameCatania) :
  Col = [[] for iCol in range(2)]
  with open(FileNameCatania) as FileCatania:
    for iLine in FileCatania:
      if "#" not in iLine:
        iLine.strip().split(',')
        iLineSep = iLine.split()
        if len(iLineSep) > 1:
          for iCol in range(2) :
            Col[iCol].append(float(iLineSep[iCol]))
  Cataniadic = {'PtCent' : Col[0],'Cent' : Col[1]}
  return Cataniadic
