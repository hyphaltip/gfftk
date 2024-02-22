from .genbank import tbl2dict, dict2tbl

def truncate(args):
    Genes = tbl2dict(args.tbl, args.fasta, debug=args.debug)
    locations = {}
    if args.cut:
        loc = args.cut.split(":")
        if loc[0] not in locations:
            locations[loc[0]] = []
        (start, end) = loc[1].split("-")
        (start, end) = (int(start), int(end))
        if start > end:
            (start, end) = (end, start)

        locations[loc[0]].append([start,end])
    # elif args.cutlist:
    #     with open(args.cutlist, "r") as f:
    #         for line in f:
    #             loc = line.strip().split(":")
    #             if loc[0] not in locations:
    #                 locations[loc[0]] = []
    #             (start, end) = loc[1].split("-")
    #             (start, end) = (int(start), int(end))
    #             if start > end:
    #                 (start, end) = (end, start)
    #             locations[loc[0]].append([start,end])
    print(locations)
    if args.cut:
        for entry in Genes:
            keep_genes = {}
            for g in entry:
                gene = entry[g]
                print(f'before gene: {gene}')
                contig = gene["contig"]
                keep_gene = 1
                if contig in locations:
                    for loc in locations[contig]:
                        print(f'loc: {loc} and gene location: {gene["location"]}')
                        print(int(gene['location'][0]))
                        print(loc)
                        if int(gene['location'][0]) >= loc[0] and int(gene['location'][0]) <= int(loc[1]):
                            # drop this feature completely it is in the truncation zone
                            print(f"deleted gene {g}")
                            keep_gene = 0
                            continue
                        elif gene['location'][0] > loc[0]:
                            # now we need to adjust the feature positions
                            
                            gene['location'] = (gene['location'][0] - loc[1] - 1, gene['location'][1] - loc[1] - 1)
                            for featuretype in ['5UTR','3UTR','CDS','mRNA']:
                                for feat_i in range(len(gene[featuretype])):
                                    gene[featuretype][feat_i] = (gene[featuretype][feat_i][0] - loc[1] - 1, gene[featuretype][feat_i][1] - loc[1] - 1)
                print(f'after gene: {gene}')
                if keep_gene:
                    keep_genes[g] = gene
            entry = keep_genes
    #dict2tbl(Genes, output=f'{args.out}.tbl')
