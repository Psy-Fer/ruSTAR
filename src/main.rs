use clap::Parser;

use ruSTAR::params::Parameters;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let params = Parameters::parse();
    ruSTAR::run(&params)
}
