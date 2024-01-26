#!/usr/bin/env julia

# Julia script to extract sequences from a fasta file for GTF features
#
# Source code: https://github.com/PoisonAlien/gtf2fasta
#
# MIT License
# Copyright (c) 2017 Anand Mayakonda <anandmt3@gmail.com>

VERSION="0.1.0"

usage = "
gtf2fasta : Extract sequences from a fasta file for GTF features
Version   : $VERSION

Usage:
gtf2fasta <input.gtf> <input.fa>
"

#--------------------- Structures
mutable struct faObject
  fa::IO
  faChr::String
  faBuf::Int
  faBufLocus::Int64
end

mutable struct bedRec
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

mutable struct Transcript
  chr::String
  txStrand::String
  geneId::String
  geneName::String
  txId::String
  exonStart::Vector{Int64}
  exonEnd::Vector{Int64}
end

function Transcript(chr::String, txStrand::String, geneId::String, geneName::String, txId::String, exonStart::Vector{Int64}, exonEnd::Vector{Int64})
    return Transcript(chr, txStrand, geneId, geneName, txId, exonStart, exonEnd)
end

#--------------------- Functions
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
        @warn("Index file not found: $fai Trying to create one with samtools.")
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
  
  if !haskey(fai, bed.bedChr)
      println(stderr, "Chromosome ", bed.bedChr, " does not exist!")
      return
  end  
  
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
        @warn("Record $chr:$startpos-$endpos exceeds chromosome length $chromLen. Skipping..")
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

    if !haskey(fai, ts.chr)
      println(stderr, "Chromosome ",ts.chr, " for ", ts.txId, " does not exist in the fasta file. Skipping..")
      return
    end  

  txFa::String = ""
  if length(ts.exonStart) == 0
    @warn("Empty Transcript with no Exons: $ts. Skipping..")
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
      txFa = txFa * extractFasta(fasta, bedRec(ts.chr, ts.exonStart[i], ts.exonEnd[i]), false, false)
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

makeExonDict = function(txrec::String)
    txdict = Dict()
    txrecSpl = split(txrec, ";")
    for i in 1:length(txrecSpl)-1
        t = split(convert(String, txrecSpl[i]), " ")
        txdict[t[1]] = t[2]
    end
    txdict
end

maketxDb = function(gtf::String)
  gtf = open(gtf, "r")

  tx = Transcript("", "", "NA", "", "", [], [])
  txDb = Dict{String, Transcript}()
  
  for ln in eachline(gtf)
      if(ln[1] != '#')
        gffspl = split(chomp(ln), "\t")
        if gffspl[3] == "exon"
            gfftx = replace(replace(gffspl[9], "\"" => ""), "; " => ";")
            txdict = makeExonDict(gfftx)
            
            if length(txDb) == 0
                tx.chr = gffspl[1]
                tx.txStrand = gffspl[7]
                tx.geneId = txdict["gene_id"]
                tx.txId = txdict["transcript_id"]
                if haskey(txdict, "gene_name")
                  tx.geneName = txdict["gene_name"]
                end
                push!(tx.exonStart, parse(Int64, gffspl[4]))
                push!(tx.exonEnd, parse(Int64, gffspl[5]))
                txDb[tx.txId] = tx
            elseif tx.txId != txdict["transcript_id"]
                if haskey(txDb, txdict["transcript_id"])
                    #If GTF is not sorted or has duplicates, throw a warning
                    @warn txdict["transcript_id"], " occurs again! Is it a sorted GTF?"
                end
                txDb[tx.txId] = tx
                tx = Transcript("", "", "NA", "", "", [], [])
                tx.chr = gffspl[1]
                tx.txStrand = gffspl[7]
                tx.geneId = txdict["gene_id"]
                tx.txId = txdict["transcript_id"]
                if haskey(txdict, "gene_name")
                  tx.geneName = txdict["gene_name"]
                end
                push!(tx.exonStart, parse(Int64, gffspl[4]))
                push!(tx.exonEnd, parse(Int64, gffspl[5]))
                txDb[tx.txId] = tx
            else
                push!(tx.exonStart, parse(Int64, gffspl[4]))
                push!(tx.exonEnd, parse(Int64, gffspl[5]))
            end
        end
      end
  end

  println(stderr, "Extracted ", length(txDb), " transcripts")
  return(txDb)
end

#--------------------- Main

if(length(ARGS) < 2)
  println(usage)
  exit()
end

gtf_file = ARGS[1]
fa_file = ARGS[2]

if !isfile(gtf_file)
  println(stderr, "GTF file $gtf_file does not exists!")
  exit()
end

if !isfile(fa_file)
  println(stderr, "Fasta file $fa_file does not exists!")
  exit()
end

gtf_basename = basename(gtf_file) * ".transcript.dict.tsv"

#Open fasta file and create an faObject
const fai = readFai(fa_file)
fa = open(fa_file, "r");
fasta = faObject(fa, "c1", 0, 0)

#Read GTF and make trasncript database
txdb = maketxDb(gtf_file)

#Write transcript database to a tsv file
open(gtf_basename, "w") do f
  hdr = "transcript_id\tgene_id\tgene_name\tstrand\tn_exons\n"
  write(f, hdr)
  for (tx, db) in txdb
    txline = db.txId * "\t" * db.geneId * "\t" * db.geneName * "\t" * db.txStrand * "\t" * string(length(db.exonStart)) * "\n"
    write(f, txline)
  end
end

#Extract fasta sequeunce for every transcript and print to console
for (tx, db) in txdb
  extractTxFasta(db, fasta)
end

println(stderr, "Wrote transcript database to: ", gtf_basename)