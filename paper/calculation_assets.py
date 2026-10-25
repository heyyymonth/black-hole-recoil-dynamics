"""Typeset recorded host and collision calculations; no physical execution."""
from pathlib import Path
import argparse,csv,hashlib,json
ROOT=Path(__file__).resolve().parents[1];PAPER=ROOT/'paper'
INPUTS=['paper/calculation_assets.py','paper/calculation_figures.py','analysis/host_reconstruction/requirements-lock.txt','results/calculations.json','results/host-reconstruction/summary.json','results/host-reconstruction/host.npz','results/host-reconstruction/species.csv','results/host-reconstruction/published-profile-comparison.csv','results/local-transport-kernels/summary.json','results/local-transport-kernels/states.csv','results/local-transport-kernels/species.csv']
OUTPUTS=['tables/host-validation.tex','tables/local-kernels.tex','figures/host-projections.pdf','figures/scattering-composition.pdf']

def rows(name):
    with (ROOT/name).open(newline='') as f:return list(csv.DictReader(f))

def tables():
    host=json.loads((ROOT/'results/host-reconstruction/summary.json').read_text())
    kernel=json.loads((ROOT/'results/local-transport-kernels/summary.json').read_text())
    def scientific(value):
        mantissa,exponent=f'{value:.2e}'.split('e')
        return '$'+mantissa+r'\times10^{'+str(int(exponent))+'}$'
    checks=host['checks']
    labels={'integrated_mass_relative_error':'Integrated mass','half_mass_radius_relative_error':'Half-mass radius','force_max_relative_error':'Force versus enclosed mass','potential_max_error_over_central_relative_potential':'Relative potential','virial_relative_error':'Virial balance','DF_max_relative_error':'DF density and raw moments','projection_max_relative_difference':'Projected luminous tracer'}
    body='\\begin{table}[t]\n\\caption{Numerical checks for the fixed NGC 3201 reconstruction. Errors are relative, except the potential error is normalized to its central value. Force checks cover $10^{-3}<r/r_h<10$; DF checks use four radii. Projected quantities refer to the luminous tracer, not empty remnant tails. These checks establish numerical consistency, not physical calibration.}\\label{tab:host-validation}\n\\centering\\begin{tabular}{lr}\n\\toprule\nCheck & Maximum error \\\\\n\\midrule\n'
    for key,label in labels.items():body+=label+' & '+scientific(checks[key])+r' \\'+'\n'
    body+='Drift quadrature refinement & '+scientific(kernel['max_A_refinement_error'])+r' \\'+'\nDiffusion quadrature refinement & '+scientific(kernel['max_D_refinement_relative_error'])+r' \\'+'\n\\bottomrule\n\\end{tabular}\n\\end{table}\n'
    states=rows('results/local-transport-kernels/states.csv');species=rows('results/local-transport-kernels/species.csv')
    body2=r'\begin{table}[t]'+'\n'+r'\caption{Local weak-scattering coefficients for a $60\,M_\odot$ diagnostic remnant with radial speed $0.25\sqrt{2\Psi(r)}$. Both drift magnitude $a_\parallel=-\mathbf A\cdot\widehat{\mathbf v}$ and diffusion trace are divided by a common constant $\ln\Lambda$. Units are $(\mathrm{km\,s^{-1}})/\mathrm{Myr}$ and $(\mathrm{km\,s^{-1}})^2/\mathrm{Myr}$, respectively. The BH fraction is its contribution to the diffusion trace in the same host. No physical impact-parameter cutoff is selected.}\label{tab:local-kernels}'+'\n'+r'\centering\begin{tabular}{rrrr}'+'\n'+r'\toprule $r/r_h$ & $a_\parallel/\ln\Lambda$ & $\mathrm{tr}\mathsf D/\ln\Lambda$ & BH fraction \\'+'\n'+r'\midrule'+'\n'
    for r in states:
        if float(r['v_over_cutoff'])!=.25 or r['orientation']!='radial':continue
        fraction=sum(float(x['trace_D_fraction']) for x in species if x['state']==r['state'] and x['type']=='BH')
        body2+=f"{float(r['r_over_rh']):g} & "+' & '.join(scientific(v) for v in [float(r['drag_over_lnLambda_kms_Myr']),float(r['trace_D_over_lnLambda_kms2_Myr']),fraction])+r' \\'+'\n'
    body2+='\\bottomrule\n\\end{tabular}\n\\end{table}\n'
    return {OUTPUTS[0]:body,OUTPUTS[1]:body2}

def checksums(paths,base):return {p:hashlib.sha256((base/p).read_bytes()).hexdigest() for p in paths}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    manifest=PAPER/'calculation-assets.json';rendered=tables()
    if args.verify:
        data=json.loads(manifest.read_text())
        assert data['input_sha256']==checksums(INPUTS,ROOT),'Calculation source changed'
        assert data['output_sha256']==checksums(OUTPUTS,PAPER),'Calculation asset changed'
        assert all((PAPER/k).read_text()==v for k,v in rendered.items()),'Table differs from evidence'
        print('Verified two calculation tables and two vector figures.');return
    for path,content in rendered.items():(PAPER/path).write_text(content)
    from calculation_figures import render
    render(ROOT,PAPER/'figures')
    manifest.write_text(json.dumps({'scope':'fixed host reconstruction and local coefficient diagnostics','input_sha256':checksums(INPUTS,ROOT),'output_sha256':checksums(OUTPUTS,PAPER)},indent=2)+'\n')
if __name__=='__main__':main()
