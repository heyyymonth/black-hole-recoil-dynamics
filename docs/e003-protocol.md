# E003 protocol: density exposure versus close-passage optical depth

**Researcher:** Venkata Hemanth Athota  
**Protocol recorded:** 19 October 2026, before execution of E003.

## Question and connection to the research problem

Does the normalized density exposure from E002 approximate the normalized opportunity for a recoiling remnant to pass near a specified population of binaries once relative velocities and gravitational focusing are included?

The larger assembly question requires retention, companion acquisition, binary evolution, and another merger. This experiment tests the connection between the retained orbit and access to potential companions. It does not simulate acquisition, exchange, hardening, or merger outcomes. It tests an approximation, not publication-level novelty.

## Declared model and comparison

The [configuration](../experiments/e003.json) is fixed before inspecting E003 outputs. Use the E001/E002 static Plummer potential and ballistic central kicks. Set the host mass to one million solar masses and its scale radius to 1 pc. The remnant mass is 60 solar masses. A dynamically negligible tracer population contains 100 binaries of total mass 60 solar masses each, with semimajor axis 0.1 AU. Their number density follows the untruncated Plummer profile and their center-of-mass velocities follow its isotropic equilibrium distribution. This is an explicitly unsegregated idealization, not a calibrated black-hole subsystem.

Use a stationary central control and kick fractions 0.3, 0.8, and 0.9. All have the same horizon, ten host dynamical times (approximately 0.149 Myr). Admit only test orbits with apocenters within five Plummer half-mass radii. This is a screening convention inspired by [Islam, section V.1.3](https://arxiv.org/html/2604.04546v2), not a self-consistently truncated host. The former 0.99 near-escape example does not pass this screen.

Define a close passage as a two-body monopole pericenter at or inside 1 AU of a target binary's center of mass. The radius is ten binary semimajor axes: a gateway to the binary's neighborhood, not a cross section for exchange or a strong three-body outcome. The isolated two-body approximation is local; it neglects host tides during the close approach and the binary's internal multipoles.

For local relative speed w, the cross section is


$$
\Sigma(w)=\pi R^2[1+2G(m_{\rm rem}+m_{\rm bin})/(R w^2)].
$$


Compute the optical depth along the prescribed orbit,


$$
\tau(T)=\int_0^T n_{\rm bin}[r(t)]\langle\Sigma(w)w\rangle\,dt.
$$


This is the kinetic-theory expected passage count in a fixed dilute bath along an unperturbed path. We do not convert it into a merger probability or assume the orbit remains unchanged after an actual encounter. Retain small absolute counts honestly; do not tune the tracer abundance to produce events.

Compare D = I(H)/H with Q = tau(eta,H)/tau(0,H). The signed proxy error is D/Q - 1, using the velocity-weighted result as denominator. Also retain the geometric and focusing contributions separately to explain any discrepancy.

## Decision rule fixed before results

- A nonzero-kick case fails this benchmark's screening tolerance if |D/Q - 1| > 0.10.
- This 10% tolerance is an operational approximation-test threshold. It is not a statistical significance level or a requirement derived from a GW inference model.
- If all cases pass, report adequacy only for these cases and stop extending this grid.
- If any case fails, identify the velocity term responsible, quantify the discrepancy, and stop this benchmark after numerical validation. Do not launch an unrestricted parameter sweep or claim failure throughout all clusters.
- Use independent angular/speed integration, a separately integrated trajectory, stationary/fast-motion analytical limits, and tighter quadrature. Require agreement better than 1e-6 relative for the reported nonzero integrals, comfortably below the screening tolerance.

## Relevance and limitations checked before execution

The exact Plummer distribution has f(E) proportional to (-E)^(7/2); see the [versioned galpy distribution documentation](https://docs.galpy.org/en/v1.12.0/reference/dfplummer.html). With y = u/v_escape(r), its local speed density is 512/(7 pi) y^2 (1-y^2)^(7/2), 0 <= y <= 1. The implementation derives its velocity moments rather than importing galpy.

[Leigh et al.](https://arxiv.org/abs/1711.10494) distinguish dynamical encounters and subsequent binary hardening toward merger. E003 isolates a simpler passage-rate ingredient. [Islam](https://arxiv.org/abs/2604.04546v2) already considers remnant displacement and recentering; E003 neither implements nor falsifies its frictional repeat-merger prescription.

The static host, unsegregated tracer bath, and unperturbed orbit are assumptions. The short horizon does not establish that dynamical friction, diffusion, collective core response, or host evolution are negligible in a real cluster. There is no observed GW dataset, stochastic realization, or measured cluster calibration here. More kick samples cannot remedy these physical omissions.
