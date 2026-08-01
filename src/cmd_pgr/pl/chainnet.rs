use clap::{ArgMatches, Command};
use cmd_lib::run_cmd;

/// Build the clap subcommand for chainnet.
pub fn make_subcommand() -> Command {
    Command::new("chainnet")
        .about("Runs the native pgr chain-net-axt-maf pipeline")
        .after_help(
            r###"
This command runs the pairwise genome alignment pipeline (psl-chain-net-axt-maf)
entirely with native pgr commands — no external kent-tools required.

* Output has been verified byte-for-byte identical against the UCSC kent-tools
  pipeline (`pgr pl ucsc`) for all intermediate files (chain, net, axt).

* <target> and <query> are fasta files
* <psl> can be a .psl file or a directory containing multiple .psl files
* Default names of target and query in the output .maf are derived from the
  basename of <target> and <query>

* `--gap-model` and `--min-score`:
    * Human18vsChimp2 use `loose` and 1000
    * Human19vsChimp3 use `medium` and 5000
    * `loose` corresponds to chicken/human linear gap costs
    * `medium` corresponds to mouse/human linear gap costs

* `--syn`: generate syntenic alignments (netFilter --syn + chainSplit path)

Pipeline steps (all native pgr):

1. Prep:    `pgr fa size` + `pgr fa to-2bit`
2. Chain:   `pgr psl chain` + `pgr chain anti-repeat`
3. Merge:   `pgr chain sort`
4. PreNet:  `pgr chain pre-net`
5. Net:     `pgr chain net` + `pgr net syntenic` + `pgr net subset` + `pgr chain stitch` + `pgr net split`
6. Axt:     `pgr net to-axt` | `pgr axt sort`
7. Maf:     `pgr axt to-maf`

Definitions:

* The *target* is the reference genome sequence
* The *query* is some other genome sequence

* A *chain* is a sequence of non-overlapping gapless blocks, with single- or
  double-sided gaps between blocks. Within a chain, target and query coords are
  monotonically non-decreasing.
* A *net* is a hierarchical collection of chains.

References:

* [Chains Nets](https://genomewiki.ucsc.edu/index.php/Chains_Nets)

"###,
        )
        .arg(crate::cmd_pgr::args::target_genome_arg(
            "Path to the target genome FA file",
        ))
        .arg(crate::cmd_pgr::args::query_genome_arg(
            "Path to the query genome FA file",
        ))
        .arg(crate::cmd_pgr::args::psl_positional_arg(
            "Path to the PSL file or directory containing PSL files",
        ))
        .arg(crate::cmd_pgr::args::gap_model_arg(
            "loose",
            &["loose", "medium"],
            "Linear gap cost setting for psl chain",
        ))
        .arg(crate::cmd_pgr::args::min_score_arg("1000"))
        .arg(crate::cmd_pgr::args::t_name_arg(None))
        .arg(crate::cmd_pgr::args::q_name_arg(None))
        .arg(crate::cmd_pgr::args::syn_arg("Generate syntenic alignments"))
        .arg(crate::cmd_pgr::args::outdir_arg())
}

/// Execute the chainnet command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outdir = args.get_one::<String>("outdir").unwrap();
    if outdir != "stdout" {
        std::fs::create_dir_all(outdir)?;
    }

    let opt_gap_model = args.get_one::<String>("gap_model").unwrap();
    let opt_minscore = *args.get_one::<f64>("min_score").unwrap();

    let is_syn = args.get_flag("syn");

    let ctx = pgr::libs::pl::PipelineCtx::new("pgr_chainnet_")?;
    let pgr = ctx.pgr.clone();

    run_cmd!(info "==> Absolute paths")?;
    let abs_target = ctx.abs_path(args.get_one::<String>("target").unwrap())?;
    let abs_query = ctx.abs_path(args.get_one::<String>("query").unwrap())?;

    let opt_tname = if let Some(tname) = args.get_one::<String>("t_name") {
        if tname.is_empty() {
            "".to_string()
        } else {
            format!("{}.", tname)
        }
    } else {
        format!("{}.", pgr::libs::io::basename_or_err(&abs_target)?)
    };
    let opt_qname = if let Some(qname) = args.get_one::<String>("q_name") {
        if qname.is_empty() {
            "".to_string()
        } else {
            format!("{}.", qname)
        }
    } else {
        format!("{}.", pgr::libs::io::basename_or_err(&abs_query)?)
    };

    let abs_psl = ctx.abs_path(args.get_one::<String>("psl").unwrap())?;
    let infiles = if std::path::Path::new(&abs_psl).is_dir() {
        pgr::libs::io::list_files_ext(&abs_psl, "psl")
    } else {
        vec![abs_psl]
    };

    let abs_outdir = pgr::libs::pl::abs_path_or_stdout(outdir)?;

    let _cwd_guard = ctx.enter()?;

    // 1. Prep: sizes and 2bit for both target and query
    run_cmd!(info "==> Target .sizes and .2bit")?;
    run_cmd!(
        ${pgr} fa size ${abs_target} -o target.chr.sizes;
        ${pgr} fa to-2bit ${abs_target} -o target.chr.2bit;
    )?;
    run_cmd!(info "==> Query .sizes and .2bit")?;
    run_cmd!(
        ${pgr} fa size ${abs_query} -o query.chr.sizes;
        ${pgr} fa to-2bit ${abs_query} -o query.chr.2bit;
    )?;

    // 2. Chain: psl chain + anti-repeat (per input PSL file)
    run_cmd!(info "==> psl chain + anti-repeat")?;
    std::fs::create_dir_all("pslChain")?;
    for infile in &infiles {
        let stem = pgr::libs::io::basename_or_err(infile)?;
        run_cmd!(
            ${pgr} psl chain target.chr.2bit query.chr.2bit ${infile}
                --min-score ${opt_minscore} --gap-model ${opt_gap_model}
                -o pslChain/${stem}.tmp
        )?;
        run_cmd!(
            ${pgr} chain anti-repeat
                --target-2bit target.chr.2bit --query-2bit query.chr.2bit
                pslChain/${stem}.tmp -o pslChain/${stem}.chain
        )?;
    }

    // 3. Merge: chain sort (all .chain files → all.chain)
    run_cmd!(info "==> chain sort")?;
    {
        let chain_files = pgr::libs::io::list_files_ext("pslChain", "chain");
        let mut tmp = tempfile::NamedTempFile::new()?;
        {
            use std::io::Write;
            for s in &chain_files {
                writeln!(tmp, "{}", s)?;
            }
            tmp.flush()?;
        }
        let tmp_path = tmp
            .path()
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("temp file path is not valid UTF-8"))?
            .to_string();
        run_cmd!(
            ${pgr} chain sort --input-list ${tmp_path} -o all.chain
        )?;
    }

    // 4. PreNet: chain pre-net
    run_cmd!(info "==> chain pre-net")?;
    run_cmd!(
        ${pgr} chain pre-net all.chain target.chr.sizes query.chr.sizes -o all.pre.chain
    )?;

    // 5. Net: chain net + net syntenic + net subset + chain stitch + net split
    run_cmd!(info "==> chain-net")?;
    {
        // chainNet: write target.net and query.net
        run_cmd!(
            ${pgr} chain net all.pre.chain target.chr.sizes query.chr.sizes
                target.chainnet query.chainnet
        )?;

        // netSyntenic: add synteny info to target net
        run_cmd!(
            ${pgr} net syntenic target.chainnet -o noClass.net
        )?;

        // netChainSubset + chainStitchId
        run_cmd!(
            ${pgr} net subset noClass.net all.chain subset.chain;
            ${pgr} chain stitch subset.chain -o over.chain;
        )?;

        // netSplit: split net into per-chromosome files
        std::fs::create_dir_all("net")?;
        run_cmd!(
            ${pgr} net split noClass.net -o net
        )?;
    }

    // 6. Axt: net to-axt | axt sort (per chromosome net file)
    run_cmd!(info "==> net to-axt + axt sort")?;
    {
        std::fs::create_dir_all("axtNet")?;
        let files = pgr::libs::io::list_files_ext("net", "net");
        for file in &files {
            let stem = pgr::libs::io::basename_or_err(file)?;
            run_cmd!(
                ${pgr} net to-axt ${file} all.pre.chain target.chr.2bit query.chr.2bit -o stdout |
                    ${pgr} axt sort stdin -o axtNet/${stem}.axt
            )?;
        }
    }

    // 7. Maf: axt to-maf (or synteny mode)
    run_cmd!(info "==> axt to-maf")?;
    if !is_syn {
        let files = pgr::libs::io::list_files_ext("axtNet", "axt");
        for file in &files {
            let stem = pgr::libs::io::basename_or_err(file)?;
            let maf_output = if abs_outdir == "stdout" {
                "stdout".to_string()
            } else {
                format!("{}/{}.maf", abs_outdir, stem)
            };
            if opt_tname.is_empty() {
                run_cmd!(
                    ${pgr} axt to-maf ${file}
                        -t target.chr.sizes -q query.chr.sizes -o ${maf_output}
                )?;
            } else {
                run_cmd!(
                    ${pgr} axt to-maf ${file}
                        --t-prefix ${opt_tname} --q-prefix ${opt_qname}
                        -t target.chr.sizes -q query.chr.sizes -o ${maf_output}
                )?;
            }
        }
    } else {
        std::fs::create_dir_all("synNet")?;

        run_cmd!(info "==> synNet.maf")?;

        // netFilter --syn + netSplit
        run_cmd!(
            ${pgr} net filter noClass.net --syn -o synNet.net;
            ${pgr} net split synNet.net -o synNet;
        )?;

        // chainSplit: split chains by target (default)
        run_cmd!(
            ${pgr} chain split all.chain -o synNet
        )?;

        // netToAxt | axtSort | axtToMaf per chromosome
        let files = pgr::libs::io::list_files_ext("synNet", "net");
        for file in &files {
            let stem = pgr::libs::io::basename_or_err(file)?;
            let net_stem = file
                .strip_suffix(".net")
                .ok_or_else(|| anyhow::anyhow!("expected .net suffix: {}", file))?;
            let chain_file = format!("{}.chain", net_stem);
            let maf_output = if abs_outdir == "stdout" {
                "stdout".to_string()
            } else {
                format!("{}/{}.maf", abs_outdir, stem)
            };
            if opt_tname.is_empty() {
                run_cmd!(
                    ${pgr} net to-axt ${file} ${chain_file} target.chr.2bit query.chr.2bit -o stdout |
                        ${pgr} axt sort stdin -o stdout |
                        ${pgr} axt to-maf stdin
                            -t target.chr.sizes -q query.chr.sizes -o ${maf_output}
                )?;
            } else {
                run_cmd!(
                    ${pgr} net to-axt ${file} ${chain_file} target.chr.2bit query.chr.2bit -o stdout |
                        ${pgr} axt sort stdin -o stdout |
                        ${pgr} axt to-maf stdin
                            --t-prefix ${opt_tname} --q-prefix ${opt_qname}
                            -t target.chr.sizes -q query.chr.sizes -o ${maf_output}
                )?;
            }
        }
    }

    // Done
    run_cmd!(info "==> Done")?;

    Ok(())
}
