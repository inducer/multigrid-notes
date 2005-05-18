<TeXmacs|1.0.5>

<style|<tuple|article|number-long-article|skript>>

<\body>
  Multigrid Methods

  Prof. Christian Wieners, SS 2005

  Korrekturen und Vorschläge an <with|font-family|tt|kloeckner@math.uni-karlsruhe.de>.

  <\table-of-contents|toc>
    <vspace*|1fn><with|font-series|bold|math-font-series|bold|Table of
    contents> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-1><vspace|0.5fn>

    1<space|2spc>Introduction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-2>

    <with|par-left|1.5fn|1.1<space|2spc>The model problem
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-3>>

    <with|par-left|1.5fn|1.2<space|2spc>Iterative Methods for Linear Systems
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-4>>

    <with|par-left|1.5fn|1.3<space|2spc>The Multigrid Idea
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-5>>

    2<space|2spc>A two-level method <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-6>

    <with|par-left|1.5fn|2.1<space|2spc>The model problem
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-7>>

    <with|par-left|1.5fn|2.2<space|2spc>Implementation of the two-level
    method <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-11>>

    3<space|2spc>Classical two-level analysis
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-12>
  </table-of-contents>

  <include|chapter01.tm>

  <include|chapter02.tm>

  <include|chapter03.tm>
</body>

<\initial>
  <\collection>
    <associate|language|american>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|<uninit>|3>>
    <associate|auto-10|<tuple|2.3|6|chapter02.tm>>
    <associate|auto-11|<tuple|2.2|10|chapter02.tm>>
    <associate|auto-12|<tuple|2.4|11|chapter02.tm>>
    <associate|auto-13|<tuple|3|11|chapter03.tm>>
    <associate|auto-2|<tuple|1|3|chapter01.tm>>
    <associate|auto-3|<tuple|1.1|3|chapter01.tm>>
    <associate|auto-4|<tuple|1.2|4|chapter01.tm>>
    <associate|auto-5|<tuple|1.3|4|chapter01.tm>>
    <associate|auto-6|<tuple|2|5|chapter02.tm>>
    <associate|auto-7|<tuple|2.1|5|chapter02.tm>>
    <associate|auto-8|<tuple|2.1|5|chapter02.tm>>
    <associate|auto-9|<tuple|2.2|6|chapter02.tm>>
    <associate|eq:twolevel-step2|<tuple|2.1|10|chapter02.tm>>
    <associate|lem:2step-energy-estimate|<tuple|2.10|8|chapter02.tm>>
    <associate|subsec:fd-model-probelm|<tuple|1.1|?>>
    <associate|subsec:fd-model-problem|<tuple|1.1|3|chapter01.tm>>
    <associate|the:twolevel|<tuple|2.9|8|chapter02.tm>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|The triangle transformation.|<pageref|auto-8>>

      <tuple|normal|A uniform grid on <with|mode|<quote|math>|(0,1)<rsup|2>>.|<pageref|auto-9>>

      <tuple|normal|Red refinement.|<pageref|auto-10>>

      <tuple|normal|Interpolation to a finer
      discretization.|<pageref|auto-12>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Table
      of contents> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|1.1<space|2spc>The Model Problem
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1.5fn>|1.2<space|2spc>Iterative Methods for
      Linear Systems <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1.5fn>|1.3<space|2spc>The Multigrid Idea
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>A
      Two-Level Method> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|2.1<space|2spc>The Model Problem
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|1.5fn>|2.2<space|2spc>Implementation of the
      Two-Level Method <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Classical
      Two-Level Analysis> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>