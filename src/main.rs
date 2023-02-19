use fancy_regex::Regex;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufReader, Read};
use std::process::exit;
use std::{env, io};
#[macro_use] // Attribute
extern crate prettytable;
use prettytable::Table;
/// Report struct contains meta data of sequence
/// In COVID19 version, default values are assigned except seqid.
#[derive(Debug)] // Attribute drive debug 主要是用于fmt方式打印，默认struct是不能的。
struct Report {
    seqid: String,
    organism: String,
    gcode: String,
    moltype: String,
    topology: String,
    strand: String,
}

impl Report {
    // implement a methods (aka. constructor) for struct
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

/* Report is printable!
similar as __str__ in python
*/
impl Display for Report {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Seqid: '{}',  Organism: '{}',  gcode: '{}',  moltype: '{}',  topology: '{}',  strand: '{}'",
            self.seqid, self.organism, self.gcode, self.moltype, self.topology, self.strand
        )
    }
}

fn validate_seq_len(
    seqid: &String,
    seq_len: &i32, /*借用本质上也是引用，所以需要用*解引用获得值 */
    min_len: &i32,
    max_len: &i32,
    msg_table: &mut Table,
) -> bool {
    if *seq_len == -1i32 {
        /* == 没有实现自动解引用，所以就没办法直接比较，必须要解引用才可以。 */
        return true;
    }
    if seq_len < min_len || seq_len > max_len {
        msg_table.add_row(row!["Sequence",format!( "For SARS-CoV-2 submission, sequence length must be between {min_len} and {max_len}. SeqLength of '{seqid}' is {seq_len}")]);
        return false;
    } else {
        return true;
    }
}

fn validate_seq_id(seqid: &String, msg_table: &mut Table) -> bool {
    let mut flag = true;
    match seqid.chars().nth(1) {
        Some(first_char) => {
            if !first_char.is_alphabetic() {
                msg_table.add_row(row!["Defline",format!( "Found invalid character '{first_char}' in seqid '{seqid}'. Seqid must starts with a letter.")]);
                flag = false;
            }
            for c in seqid[1..].chars() {
                if !(c.is_numeric()
                    || c.is_alphabetic()
                    || c == '_'
                    || c == '*'
                    || c == '#'
                    || c == '.'
                    || c == '-'
                    || c == ':')
                {
                    msg_table.add_row(row!["Defline",format!( "Found invalid character: '{c}' in seqid: '{seqid}'. Only letters, numbers, '_', '-', '*', '#', '.', ':' are permitted.")]);
                    flag = false;
                }
            }
            let seqid_len = seqid.len();
            if seqid_len > 24 {
                msg_table.add_row(row![
                    "Defline",
                    format!(
                        "Seqid max length is 23, found length of {seqid_len} for seqid: {seqid}."
                    )
                ]);
                flag = false;
            }
        }
        None => {
            msg_table.add_row(row![
                "Defline",
                format!("Found invalid seqid '{seqid}'. seqid must has at least one character.")
            ]);
            return false;
        }
    }
    return flag;
}

fn validate_seq_n_pct(
    seq_len: &i32,
    previous_seqid: &String,
    seq_n_count: &i32,
    msg_table: &mut Table,
) {
    if *seq_len > 0 {
        // Rust有自动解引用的方法的，但是有的时候需要手动。
        let seq_n_pct = (*seq_n_count as f32) / (*seq_len as f32);
        if seq_n_pct >= 0.5 {
            let seq_n_pct_scaled = seq_n_pct * 100.0;
            msg_table.add_row(row!["Defline",format!( "For SARS-CoV-2 submission, the proportion of unknown bases in the sequence exceeds 50% is not allowed. Found {seq_n_count}/{seq_len}({seq_n_pct_scaled}%) for sequence '{previous_seqid}' ")]);
        }
    }
}

fn validate_seqid_unique(report_list: &Vec<Report>, msg_table: &mut Table) {
    let mut seqid_map: HashMap<&str, &str> = HashMap::new();

    for x in report_list.iter() {
        let seqid_value = x.seqid.as_str();
        if seqid_map.contains_key(seqid_value) {
            msg_table.add_row(row![
                "Defline",
                format!("Found duplicated sequence id: '{seqid_value}'")
            ]);
        } else {
            seqid_map.insert(seqid_value, "");
        }
    }
}

/// Get sequence id from defline
/// Using fancy_regex to enable look-head and look-behind
fn get_seqid(defline: &String) -> String {
    let re = Regex::new(r"^(>.*?)[\s|\n]").unwrap();
    let caps = re.captures(defline);
    match caps {
        Ok(t) => match t {
            Some(x) => {
                return x.get(1).unwrap().as_str().to_string();
            }
            None => {
                return "".to_string();
            }
        },
        Err(_e) => {
            return "".to_string();
        }
    }
}

/*
    实现动态类型的思路
    1. dyn trait 这个trait object的升级，是运行时的动态分发。跟泛型不一样，泛型是编译器展开，是静态多态。
    2. 实现函数返回值的多个类型
        通过enum来包装多个类型
        函数返回值是Result<enumK,&str>
        函数外部 使用 match enum类型变量，同时使用 let v = match enumK {} 方式来返回不同类型。
*/

const BUF_SIZE: usize = 1024 * 1024 * 4;
trait Readbuf {
    fn read_buff(&mut self, buf: &mut [u8; BUF_SIZE]) -> io::Result<usize>;
}

struct Gzreader {
    reader: BufReader<flate2::read::GzDecoder<File>>,
}

impl Readbuf for Gzreader {
    fn read_buff(&mut self, buf: &mut [u8; BUF_SIZE]) -> io::Result<usize> {
        self.reader.read(buf)
    }
}

struct Txtreader {
    reader: File,
}

impl Readbuf for Txtreader {
    fn read_buff(&mut self, buf: &mut [u8; BUF_SIZE]) -> io::Result<usize> {
        self.reader.read(buf)
    }
}

fn get_fa_reader2(_fa_file: &String) -> Box<dyn Readbuf> {
    /*函数最好返回一个Result */
    let file = File::open(_fa_file).expect("Some Error Occurred, can not open Fasta file.");
    let binding = _fa_file.split(".").last().unwrap().to_lowercase();
    let suffix = binding.as_str();
    if suffix == "gz" {
        let gz = Gzreader {
            reader: BufReader::new(GzDecoder::new(file)),
        };
        return Box::new(gz);
    } else if (suffix == "fa") || (suffix == "fsa") || (suffix == "fna") || (suffix == "fasta") {
        let txt = Txtreader { reader: file };
        return Box::new(txt);
    } else {
        /*必须有这一个分支处理其他情况否则报错 */
        eprintln!("The suffix of sequence file should be one of [.fa, .fsa, .fna, .fasta] or it should be a combination of those with .gz ");
        exit(1)
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    //Validate the file path
    if args.len() < 2 {
        println!("Please specify Fasta file");
        return;
    }

    /* use pretty table to store messages */
    let mut msg_table = Table::new();
    // msg_table.add_row(row!["Error_type", "Message"]);
    msg_table.set_titles(row!["Error_type", "Message"]);

    /* define constant values */
    const MIN_LEN: i32 = 50;
    const MAX_LEN: i32 = 30_000;

    let _fa_file: &String = &args[1];
    let mut buf = [0u8; BUF_SIZE];
    /* compatible with flat file or gz format */
    /* 这里有必要存在这个变量，因为只有执行了这一句，把get_fa_reader返回值了，
    才能一次性打开文件，如果在loop中，则是会多次打开。
    */

    // 这个是适配get_fa_reader 函数的，因为函数返回值是一个Result，所以就需要解包
    // let mut file_reader = match get_fa_reader(&_fa_file) {
    //     Ok(t) => t,
    //     Err(e) => {
    //         eprintln!("{}", e);
    //         exit(1)
    //     }
    // }; /*这里直接unwrap解包 */
    /*初始化fasta reader */
    let mut file_reader2 = get_fa_reader2(&_fa_file);
    let mut previous_char: char = '@';
    let mut previous_2nd_char: char = '#';
    let mut report_list: Vec<Report> = Vec::new(); //collections type
    let mut is_in_defline: bool = false;
    let mut defline: String = String::from("");
    let mut previous_defline: String = String::from("");
    // let mut previous_seqid: String = String::from(""); // removed here to void unexpected manner when using it at iterations
    let mut line_number: i32 = 1;
    let mut seq_len: i32 = -1;
    let mut seq_n_count: i32 = 0;
    let mut chr_pos: i32 = 0;

    /*
        TODO:
        pre-define u8 code, can be used to comparison with nucleotide in u8 format directly.
    */
    // const CHAR_END:u8 = '\0' as u8;

    loop {
        /* let v = match u ... 这种模式放置了sub-scope中的内容无法再outer-scope中访问 */
        let num_bytes = file_reader2.read_buff(&mut buf);

        if num_bytes.unwrap() == 0 {
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
            /*
            初始化的buf为0u8, 因此当*c的值为 0u8的时候，
            说明这个buf没有占满，但是不需要处理了。
            当判断 c 为 0u8 也就是字符 '\0'的时候，退出循环。
            */
            if *c == ('\0' as u8) {
                break;
            }
            let cc = *c as char;
            /* 统计行数并且重置chr的位置 */
            if cc == '\n' {
                line_number += 1;
                chr_pos = 0; //换行重置为0
            }
            /* The major brach is to distinguish defline related chars from sequence related chars. */
            if cc == '>' {
                if previous_char == '@' || previous_char == '\n' {
                    is_in_defline = true;
                    defline.push(cc);
                    /* 判断序列是否为N结尾,程序靠后部分还有一段代码 */
                    let previous_seqid = get_seqid(&previous_defline);
                    if previous_2nd_char == 'N' || previous_2nd_char == 'n' {
                        msg_table.add_row(row!["Nucleotide",format!("Found invalid 'N' at end of sequence '{}'. It should not end with 'N' or 'n'", previous_seqid)]);
                    }
                } else {
                    let previous_seqid = get_seqid(&previous_defline);
                    chr_pos += 1;
                    msg_table.add_row(row!["Nucleotide",format!("Found invalid '>' at Line {}, Column {} in sequence(seqid:'{}'). This symbol is not allowed in the sequence. Please check whether the new-line character is missing.",line_number,chr_pos, previous_seqid)]);
                }
            } else {
                /* 判断defline结束，开始处理defline */
                if previous_char == '\n' && is_in_defline == true {
                    is_in_defline = false;
                    let previous_seqid = get_seqid(&previous_defline);
                    if previous_defline.len() > 0 {
                        validate_seq_id(&previous_seqid, &mut msg_table);
                        validate_seq_len(
                            &previous_seqid,
                            &seq_len,
                            &MIN_LEN,
                            &MAX_LEN,
                            &mut msg_table,
                        );
                        let mut previous_report = Report::new();
                        previous_report.seqid = previous_seqid.clone();
                        report_list.push(previous_report);
                        /* validate N percent, 循环退出后再次调用一次，处理最后的序列 */
                        validate_seq_n_pct(&seq_len, &previous_seqid, &seq_n_count, &mut msg_table);
                    }
                    /* 重置 */
                    previous_defline = defline.clone();
                    defline = String::from("");
                    seq_len = 0;
                    seq_n_count = 0;
                }
                /* 开始存储defline */
                if is_in_defline {
                    defline.push(cc);
                } else {
                    /* 这里是校验每个字符，统计长度 */
                    chr_pos += 1;
                    if cc != '\n' {
                        seq_len = seq_len + 1; // 忽略换行符
                    }
                    /* 判断开始为N的报错 */
                    if seq_len == 1 && (cc == 'N' || cc == 'n') {
                        let previous_seqid_local = get_seqid(&previous_defline);
                        msg_table.add_row(row!["Nucleotide",format!("Found invalid 'N' at start of sequence '{}'. It should not start with 'N' or 'n'", previous_seqid_local)]);
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
                        'N' => {
                            seq_n_count += 1;
                        }
                        'n' => {
                            seq_n_count += 1;
                        }
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
                            msg_table.add_row(row![
                                "Nucleotide",
                                format!(
                                    "Found invalid char '{}' at Line {}, Column {}",
                                    cc, line_number, chr_pos
                                )
                            ]);
                        }
                    }
                }
            }
            previous_2nd_char = previous_char.clone();
            previous_char = cc.clone();
        }
        /*
        必须每次迭代之后对buff重置，否则最后一次迭代无法充满
        整个buff导致buff后边的内容是非法的但是又能被程序识别，
        造成不可预知的错误。*/
        buf.fill(0u8);
    }
    /*process the final sequence */
    let previous_seqid = get_seqid(&previous_defline);
    validate_seq_id(&previous_seqid, &mut msg_table);
    validate_seq_len(
        &previous_seqid,
        &seq_len,
        &MIN_LEN,
        &MAX_LEN,
        &mut msg_table,
    );
    /* 最后一条序列的最后一个碱基存储在previous_char中 */
    if previous_2nd_char == 'N' || previous_2nd_char == 'n' {
        msg_table.add_row(row![
            "Nucleotide",
            format!(
                "Found invalid 'N' at end of sequence '{}'. It should not end with 'N' or 'n'",
                previous_seqid
            )
        ]);
    }
    validate_seq_n_pct(&seq_len, &previous_seqid, &seq_n_count, &mut msg_table);
    let mut previous_report = Report::new();
    previous_report.seqid = previous_seqid;
    report_list.push(previous_report);
    validate_seqid_unique(&report_list, &mut msg_table);
    /*print 2seqidcheck.txt */
    println!("seqid\torganism\tgenetic_code\tmoltype\ttopology\tstrand");
    for x in report_list.iter() {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            x.seqid, x.organism, x.gcode, x.moltype, x.topology, x.strand
        );
        // println!("{}",x);
    }

    if msg_table.len() > 0 {
        eprint!("{}", msg_table);
    }
}
