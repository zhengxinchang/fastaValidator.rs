use fancy_regex::Regex;
use std::env;
use std::fs;
use std::io::Read;

/// Report struct contains meta data of sequence
/// In COVID19 version, default values are assigned except seqid.
#[derive(Debug)]
struct Report {
    seqid: String,
    organism: String,
    gcode: String,
    moltype: String,
    topology: String,
    strand: String,
}

impl Report {
    fn new() -> Self {
        Report {
            seqid: String::from(""),
            organism: String::from("Severe acute respiratory syndrome coronavirus 2"),
            gcode: String::from("1"),
            moltype: String::from("genomic RNA"),
            topology: String::from("linear"),
            strand: String::from("single"),
        }
    }
}

fn validate_seq_len(seqid: &String, seq_len: &i32, min_len: &i32, max_len: &i32) -> bool {
    if *seq_len == -1i32 {
        return true;
    }
    if (seq_len < min_len) || (seq_len > max_len) {
        // eprintln!("SeqLen {}", seq_len);

        eprintln!("For SARS-CoV-2 submission, sequence length must be between {min_len} and {max_len}. SeqLength of {seqid} is {seq_len}", );
        return false;
    } else {
        return true;
    }
}

fn validate_seq_id(seqid: &String) -> bool {
    let mut flag = true;
    match seqid.chars().nth(1) {
        Some(first_char) => {
            if !first_char.is_alphabetic() {
                eprintln!("Found invalid character '{first_char}' in seqid '{seqid}'. Seqid must starts with a letter.");
                flag = false;
            }
            for c in seqid[1..].chars() {
                if !(c.is_numeric()
                    || c.is_alphabetic()
                    || (c == '_')
                    || (c == '*')
                    || (c == '#')
                    || (c == '.')
                    || (c == '-')
                    || (c == ':'))
                {
                    eprintln!("Found invalid character: '{c}' in seqid: '{seqid}'. Only letters, numbers, '_', '-', '*', '#', '.', ':' are permitted.");
                }
            }
            let seqid_len = seqid.len();
            if seqid_len > 24 {
                eprintln!("Seqid max length is 23, found length of {seqid_len} for seqid: {seqid}.");
            }
        }
        None => {
            eprintln!("Found invalid seqid '{seqid}'. seqid must has at least one character.");
            return false;
        }
    }
    return flag;
}


fn validate_seq_n_pct(seq_len: &i32,previous_seqid: &String,seq_n_count:&i32){
    // let mut seq_n_pct:f32 = 0.0;
    if *seq_len > 0 { // Rust有自动解引用的方法的，但是有的时候需要手动。
        let seq_n_pct = (*seq_n_count as f32) /(*seq_len as f32);
        if seq_n_pct >=0.50  {
            let seq_n_pct_scaled = seq_n_pct * 100.0 ;
            eprintln!(
                "For SARS-CoV-2 submission, the proportion of unknown bases in the sequence exceeds 50% is not allowed. Found {seq_n_count}/{seq_len}({seq_n_pct_scaled}%) for sequence '{previous_seqid}' "
                
            )
        }
    }
}

/// Get squence id from defline
/// Using fancy_regex to enable look-head and look-behind
fn get_seqid(defline: &String) -> String {
    let re = Regex::new(r"^(>.*?)[\s|\n]").unwrap();
    let caps = re.captures(defline);
    match caps {
        Ok(t) => {
            match t {
                Some(x) => {
                    return x.get(1).unwrap().as_str().to_string();
                }
                None => return "".to_string(),
            }
        }
        Err(_e) => {
            // println!("Can not extract seqid from defline {defline},beacuse {e}");
            return "".to_string();
        }
    };
}

fn main() {
    let args: Vec<String> = env::args().collect();
    //Validate the file path
    if args.len() < 2 {
        println!("Please specify fasta file");
        return;
    }
    const BUF_SIZE: usize = 1024  * 512 ;
    const MIN_LEN: i32 = 50;
    const MAX_LEN: i32 = 30_000;
    let _fa_file: &String = &args[1];
    let mut buf = [0u8; BUF_SIZE];
    let mut file = fs::File::open(_fa_file).expect("Some Error Occurred, Can not open Fasta file.");
    let mut previous_char: char = '@';
    let mut previous_2nd_char: char = '#';
    let mut report_list: Vec<Report> = Vec::new();
    let mut is_in_defline: bool = false;
    let mut defline: String = String::from("");
    let mut previous_defline: String = String::from("");
    let mut previous_seqid: String = String::from("");
    let mut line_number: i32 = 1;
    let mut seq_len: i32 = -1;
    let mut seq_n_count:i32 = 0;
    let mut chr_pos: i32 = 0;

    while let Ok(num_bytes) = file.read(&mut buf) {
        // println!("{}", num_bytes);
        if num_bytes == 0 {
            break;
        }

        //  ## 运行逻辑
        //
        //  将整个文件看作一个一行序列，换行符等符号都是作为普通符号。
        //  遍历字符
        //   1. 首先判断是否是defline，如果是defline则将字符存储。否则直接进行校验。
        //   2. 判断defline的依据不能仅靠单个字符，要靠当前字符和前一个字符。如果前一个字符是 '@' 当前字符是 '>' 则是第一个defline。
        //       如果是后边的defline，则是 '\n' 和 '>' 作为defline的开始的flag。
        //       如何判断defline,结束， 通过检测'\n' 符号，如果'\n'前边不是 '>',并且 is_def_line == true 则判断为defline
        //       结束。
        //   3. defline的操作
        //     3.1 提取seqid
        //     3.2 校验seqid格式
        //
        //   4. 校验长度
        //   5. 校验结尾N
        //
        for c in buf.iter() {
            /* 初始化的buf为0u8, 因此当*c的值为 0u8的时候，说明这个buf没有占满，但是不需要处理了。
            当判断 c 为 0u8 也就是字符 '\0'的时候，退出循环 */
            if *c == 0u8 {
                break;
            }
            let cc = *c as char;
            /*统计行数并且重置chr的位置 */
            if cc == '\n' {
                line_number += 1;
                chr_pos = 0; //换行重置为0
            }
            if cc == '>' {
                if (previous_char == '@') || (previous_char == '\n') {
                    is_in_defline = true;
                    defline.push(cc);

                    /*判断序列是否为N结尾,程序靠后部分还有一段代码 */
                    // println!("previous_2nd_char:{previous_2nd_char},previous_char:{previous_char},cc:{cc}");
                    if (previous_2nd_char == 'N') || (previous_2nd_char == 'n') {
                        eprintln!(
                            "Found invalid 'N' at end of sequence '{}'. It should not end with 'N' or 'n'",
                            previous_seqid
                        )
                    }
                }else {
                    eprintln!(
                        "Found invalid '>' at sequence(seqid:'{}'). This symbol is not allowed in the sequence. Please check whether the new-line character is missing.",
                        previous_seqid
                    )
                }
            } else {
                /*判断defline结束，开始处理defline */
                if (previous_char == '\n') && (cc != '>') && (is_in_defline == true) {
                    is_in_defline = false;
                    if previous_defline.len() > 0 {
                        validate_seq_id(&previous_seqid);
                        validate_seq_len(&previous_seqid, &seq_len, &MIN_LEN, &MAX_LEN);
                        let mut previous_report = Report::new();
                        previous_report.seqid = previous_seqid.clone();
                        report_list.push(previous_report);

                        /*validate N percent, 循环退出后再次调用一次，处理最后的序列 */
                        validate_seq_n_pct(&seq_len,&previous_seqid,&seq_n_count);
                    }
                    /*重置*/
                    previous_defline = defline.clone();
                    previous_seqid = get_seqid(&previous_defline);
                    defline = String::from("");
                    seq_len = 0;
                    seq_n_count = 0;
                }
                /*开始存储defline */
                if is_in_defline {
                    defline.push(cc);
                } else {
                    /*这里是校验每个字符，统计长度 */
                    chr_pos += 1;
                    if cc != '\n'{
                        seq_len += 1; // 忽略换行符
                    }
                    
                    /* 判断开始为N的报错 */
                    if (seq_len == 1) && ((cc == 'N') || (cc == 'n')) {
                        eprintln!(
                            "Found invalid 'N' at start of sequence '{}'. It should not start with 'N' or 'n'",
                            previous_seqid
                        )
                    }
                    match cc {
                        'A' => {}
                        'T' => {}
                        'C' => {}
                        'G' => {}
                        'a' => {}
                        't' => {}
                        'c' => {}
                        'g' => {}
                        'N' => {seq_n_count +=1;}
                        'n' => {seq_n_count +=1;}
                        'R' => {}
                        'r' => {}
                        'Y' => {}
                        'y' => {}
                        'M' => {}
                        'm' => {}
                        'K' => {}
                        'k' => {}
                        'S' => {}
                        's' => {}
                        'W' => {}
                        'w' => {}
                        'H' => {}
                        'h' => {}
                        'B' => {}
                        'b' => {}
                        'V' => {}
                        'v' => {}
                        'D' => {}
                        'd' => {}
                        '\n' => {} // \n必须允许，不然每次换行都报错
                        _ => {
                            eprintln!(
                                "Found invalid char '{}' at Line {}, Column {}",
                                cc, line_number, chr_pos
                            )
                        }
                    }
                }
            }
            previous_2nd_char = previous_char;
            previous_char = cc;
            // println!("{cc}");
        }
    }
    /*process the final sequence */
    let previous_seqid = get_seqid(&previous_defline);
    validate_seq_id(&previous_seqid);
    validate_seq_len(&previous_seqid, &seq_len, &MIN_LEN, &MAX_LEN);
    /* 最后一条序列的最后一个碱基存储在previous_char中 */
    if (previous_char == 'N') || (previous_char == 'n') {
        eprintln!(
            "Found invalid 'N' at end of sequence '{}'. It should not end with 'N' or 'n'",
            previous_seqid
        )
    }
    validate_seq_n_pct(&seq_len, &previous_seqid, &seq_n_count);
    let mut previous_report = Report::new();
    previous_report.seqid = previous_seqid;
    report_list.push(previous_report);

    /*print 2seqidcheck.txt */
    println!("seqid\torganism\tgenetic_code\tmoltype\ttopology\tstrand");
    for x in report_list.iter() {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            x.seqid, x.organism, x.gcode, x.moltype, x.topology, x.strand
        );
    }
}
