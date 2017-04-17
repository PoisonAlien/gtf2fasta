#!/usr/bin/env julia

usage = "Tookit for gtf format conversion and sequence extraction

gtfsanta gtf2fasta <input.gtf> <input.fa> : extract sequence from gtf file
gtfsanta bed2fasta <input.bed> <input.fa> : extract sequence from bed file
"

if(length(ARGS) < 3)
  println(usage)
  exit()
end

#--------------------- Objects and functions
type faObject
  fa::IO
  faChr::String
  faBuf::Int
  faBufLocus::Int64
end

type bedRec
    bedChr::String
    bedStart::Int64
    bedEnd::Int64

    function bedRec(bedChr::String, bedStart::Int64, bedEnd::Int64)
        if(bedStart > bedEnd)
            println("$bedChr:$bedStart-$bedEnd")
            error("End position greater than Start position")
        end
        new(bedChr, bedStart, bedEnd)
    end
end

type BED12
  chrom::String
  chromStart::Int64
  chromEnd::Int64
  name::String
  score::Int64
  strand::String
  thickStart::Int64
  thickEnd::Int64
  itemRgb::Int64
  blockCount::Int64
  blockStart::Vector{Int64}
  blockEnds::Vector{Int64}

  # function bedRec(chrom::String, chromStart::Int64, chromEnd::Int64, name::String, score::Int64, strand::String, thickStart::Int64,
  #                 thickEnd::Int64, itemRgb::Int64, blockCount::Int64, blockStart::Vector{Int64}, blockEnds::Vector{Int64})
  #     if(chromStart > chromEnd)
  #         println("$name::$chrom:$chromStart-$chromEnd")
  #         error("End position greater than Start position")
  #     end
  #     new(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockStart, blockEnds)
  # end
end

printFasta = function(faseq::String, bed::bedRec, lwd::Int64 = 60)
  println(">", bed.bedChr, ":", bed.bedStart, "-", bed.bedEnd)
  lineTrack::Int64 = 1
  nlines = Int(floor(length(faseq) / lwd))
  for i in 1:nlines
      println(faseq[lineTrack:(lwd*i)])
      lineTrack = (lwd*i)+1
  end
  println(faseq[lineTrack:length(faseq)])
end

readFai = function(fa::String)
    fai = fa * ".fai"
    if !isfile(fai)
        warn("Index file not found: $fai Trying to create one with samtools.")
        run(`samtools faidx $fa`)
        if !isfile(fai)
          error("Could not create index file: $fai Is samtools installed ?")
        end
    end

    fai = open(fai, "r")
    faDict = Dict()

    for ln in eachline(fai)
        lnspl = convert(Array{String}, split(chomp(ln), "\t"))
        faDict[lnspl[1]] = map(x -> parse(Int64, x),lnspl[2:5])
    end
    close(fai)
    faDict
end

setBuffer = function(fasta::faObject, bed::bedRec)
    bedOffsets = fai[bed.bedChr]

    if bed.bedStart < bedOffsets[3]
        seekpos = bedOffsets[2]
        fasta.faBufLocus = 0
    elseif bed.bedStart % bedOffsets[3] > 1
        seekpos = Int(((floor(bed.bedStart/bedOffsets[3])*bedOffsets[4]))+bedOffsets[2])
        fasta.faBufLocus = floor(bed.bedStart/bedOffsets[3])*bedOffsets[3]
    else
        seekpos = Int((((floor(bed.bedStart/bedOffsets[3])-1)*bedOffsets[4]))+bedOffsets[2])
        fasta.faBufLocus = (floor(bed.bedStart/bedOffsets[3])-1)*bedOffsets[3]
    end

    seek(fasta.fa, seekpos)
    fasta.faChr = bed.bedChr
    fasta.faBuf = position(fasta.fa)
end

extractFasta = function(fasta::faObject, bed::bedRec, isBed::Bool = true, printFa::Bool = true)
    setBuffer(fasta, bed)
    startpos::Int64 = 0
    chr::String, endpos::Int64 = bed.bedChr, bed.bedEnd
    chromLen::Int64 = fai[bed.bedChr][1]

    if isBed
        startpos = bed.bedStart+1 #make 1-based index for start pos
    else
        startpos = bed.bedStart
    end

    if(endpos > chromLen)
        warn("Record $chr:$startpos-$endpos exceeds chromosome length $chromLen. Skipping..")
        return
    end

    basesRead::Int64 = fasta.faBufLocus
    faseq::String = ""
    faseqStart::Int64 = 0
    faseqEnd::Int64 = 0

    for ln in eachline(fasta.fa)
        ln = chomp(ln)
        basesRead = basesRead + length(ln)

        if basesRead >= startpos
            if faseqStart == 0
                faseqStart = basesRead - length(ln)
            end
            faseq = faseq * ln
            if basesRead >= endpos
                faseqEnd = basesRead
                faseqFrom = abs(startpos-faseqStart)
                faseqTo = (faseqFrom + (endpos - startpos))
                if printFa
                  printFasta(faseq[faseqFrom:faseqTo], bed)
                else
                  return(faseq[faseqFrom:faseqTo])
                end
                break
            end
        end
    end
end

##Transcript classes
type Transcript
  chr::String
  txStart::Int64
  txEnd::Int64
  txStrand::String
  geneId::String
  txId::String
  exonStart::Vector{Int64}
  exonEnd::Vector{Int64}
end

maketxdict = function(txrec::String)
    txdict = Dict()
    txrecSpl = split(txrec, ";")
    for i in 1:length(txrecSpl)-1
        t = split(convert(String, txrecSpl[i]), " ")
        txdict[t[1]] = t[2]
    end
    txdict
end

reverseComplement = function(fa::String)
    dnaDict = Dict('A' => "T", 'T' => "A", 'G' => "C", 'C' => "G", 'N' => "N")
    fa = reverse(fa)
    faCompl::String = ""
    fal::Int64 = 1
    while fal <= length(fa)
        faCompl = faCompl * dnaDict[fa[fal]]
        fal = fal + 1
    end
    faCompl
end

extractTxFasta = function(ts::Transcript, fasta::faObject)
  txFa::String = ""
  if length(ts.exonStart) == 0
    warn("Empty Transcript with no Exons: $ts. Skipping..")
    return
  end

  #Sort exons based on start position
  tsExonOrder = sortperm(ts.exonEnd)
  ts.exonStart = ts.exonStart[tsExonOrder]
  ts.exonEnd = ts.exonEnd[tsExonOrder]

  if ts.txStrand == "-"
    for i in length(ts.exonStart):-1:1
      txFa = txFa * reverseComplement(extractFasta(fasta, bedRec(ts.chr, ts.exonStart[i], ts.exonEnd[i]), false, false))
    end
  else
    for i in 1:length(ts.exonStart)
      txFa = txFa * extractFasta(fasta, bedRec(ts.chr, tx.exonStart[i], ts.exonEnd[i]), false, false)
    end
  end

  printTxFasta(txFa, ts)
end

printTxFasta = function(faseq::String, ts::Transcript, lwd::Int64 = 70)
  println(">", ts.txId)
  lineTrack::Int64 = 1
  nlines = Int(floor(length(faseq) / lwd))
  for i in 1:nlines
      println(faseq[lineTrack:(lwd*i)])
      lineTrack = (lwd*i)+1
  end
  if lineTrack <= length(faseq)
    println(faseq[lineTrack:length(faseq)])
  end
end

gffLineParse = function(gff::String, fasta::faObject)
  gffspl = split(chomp(gff), "\t")
  if(gffspl[3] == "transcript")
    gfftx = replace(replace(gffspl[9], "\"", ""), "; ", ";")
    txdict = maketxdict(gfftx)
    if firstEntry
      global firstEntry = false
      tx.chr = convert(String, gffspl[1])
      tx.txStart = parse(Int64, gffspl[4])
      tx.txEnd = parse(Int64, gffspl[5])
      tx.txStrand = convert(String, gffspl[7])
      tx.geneId = txdict["gene_id"]
      tx.txId = txdict["transcript_id"]
    else
      extractTxFasta(tx, fasta)
      tx.chr = convert(String, gffspl[1])
      tx.txStart = parse(Int64, gffspl[4])
      tx.txEnd = parse(Int64, gffspl[5])
      tx.txStrand = convert(String, gffspl[7])
      tx.geneId = txdict["gene_id"]
      tx.txId = txdict["transcript_id"]
      tx.exonStart = []
      tx.exonEnd = []
    end
  elseif(gffspl[3] == "exon")
    tx.exonStart = push!(tx.exonStart, parse(Int64, gffspl[4]))
    tx.exonEnd = push!(tx.exonEnd, parse(Int64, gffspl[5]))
  end
end

gffExonLineParse = function(gff::String, fasta::faObject)
  gffspl = split(chomp(gff), "\t")
    if gffspl[3] == "exon"
        gfftx = replace(replace(gffspl[9], "\"", ""), "; ", ";")
        txdict = maketxdict(gfftx)
        if firstEntry
          global firstEntry = false
          tx.chr = convert(String, gffspl[1])
          #tx.txStart = parse(Int64, gffspl[4])
          #tx.txEnd = parse(Int64, gffspl[5])
          tx.txStrand = convert(String, gffspl[7])
          tx.geneId = txdict["gene_id"]
          tx.txId = txdict["transcript_id"]
          tx.exonStart = push!(tx.exonStart, parse(Int64, gffspl[4]))
          tx.exonEnd = push!(tx.exonEnd, parse(Int64, gffspl[5]))
        else
            if tx.txId != txdict["transcript_id"]
                extractTxFasta(tx, fasta)
                tx.chr = convert(String, gffspl[1])
                #tx.txStart = parse(Int64, gffspl[4])
                #tx.txEnd = parse(Int64, gffspl[5])
                tx.txStrand = convert(String, gffspl[7])
                tx.geneId = txdict["gene_id"]
                tx.txId = txdict["transcript_id"]
                tx.exonStart = []
                tx.exonEnd = []
                tx.exonStart = push!(tx.exonStart, parse(Int64, gffspl[4]))
                tx.exonEnd = push!(tx.exonEnd, parse(Int64, gffspl[5]))
            else
              tx.exonStart = push!(tx.exonStart, parse(Int64, gffspl[4]))
              tx.exonEnd = push!(tx.exonEnd, parse(Int64, gffspl[5]))
            end
        end
    end
end

makeBED12 = function(gff::String)
  gffspl = split(chomp(gff), "\t")
  if(gffspl[3] == "transcript")
    gfftx = replace(replace(gffspl[9], "\"", ""), "; ", ";")
    txdict = maketxdict(gfftx)
    if firstEntry
      global firstEntry = false
      tx.chr = convert(String, gffspl[1])
      tx.txStart = parse(Int64, gffspl[4])
      tx.txEnd = parse(Int64, gffspl[5])
      tx.txStrand = convert(String, gffspl[7])
      tx.geneId = txdict["gene_id"]
      tx.txId = txdict["transcript_id"]
    else
      extractTxFasta(tx, fasta)
      tx.chr = convert(String, gffspl[1])
      tx.txStart = parse(Int64, gffspl[4])
      tx.txEnd = parse(Int64, gffspl[5])
      tx.txStrand = convert(String, gffspl[7])
      tx.geneId = txdict["gene_id"]
      tx.txId = txdict["transcript_id"]
      tx.exonStart = []
      tx.exonEnd = []
    end
  elseif(gffspl[3] == "exon")
    tx.exonStart = push!(tx.exonStart, parse(Int64, gffspl[4]))
    tx.exonEnd = push!(tx.exonEnd, parse(Int64, gffspl[5]))
  end
end

gtf2fa = function(gtf::String, fa::String)
  fa = open(fa, "r");
  gtf = open(gtf, "r")

  #Constructors
  global tx = Transcript("chr", 1, 1, "+", "g", "t", [], [])
  fasta = faObject(fa, "c1", 0, 0)
  global firstEntry = true

  for ln in eachline(gtf)
      if(ln[1] != '#')
          gffExonLineParse(ln, fasta)
      end
  end
  extractTxFasta(tx, fasta)

  close(gtf)
  close(fa)
end

bed2fa = function(inputBed::String, fa::String)
  fa = open(fa, "r");
  inputBed = open(inputBed, "r")
  fasta = faObject(fa, "c1", 0, 0)

  for bedEntry in eachline(inputBed)
    if(bedEntry[1] != '#')
      bedspl = split(chomp(bedEntry), "\t")
      bedspl = bedRec(convert(String, bedspl[1]), parse(Int64, bedspl[2]), parse(Int64, bedspl[3]))
      extractFasta(fasta, bedspl)
    end
  end

  close(fa)
  close(inputBed)
end
#--------------------- End of Objects and functions


if(ARGS[1] == "gtf2fasta")
  const fai = readFai(ARGS[3])
  gtf2fa(ARGS[2], ARGS[3])
elseif(ARGS[1] == "bed2fasta")
  const fai = readFai(ARGS[3])
  bed2fa(ARGS[2], ARGS[3])
else
  println(usage)
  exit()
end
