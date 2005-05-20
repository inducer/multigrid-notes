<TeXmacs|1.0.5>

<project|multigrid.tm>

<style|<tuple|article|skript|number-long-article>>

<\body>
  <assign|section-nr|0><section|Introduction>

  <subsection|The Model Problem>

  <label|subsec:fd-model-problem>On <with|mode|math|\<Omega\>=(0,1)<rsup|2>>,
  consider the Laplace problem:

  <\quotation>
    Find <with|mode|math|u:<wide|\<Omega\>|\<bar\>>\<rightarrow\>\<bbb-R\>>
    such that

    <\eqnarray*>
      <tformat|<table|<row|<cell|-\<Delta\>u>|<cell|=>|<cell|f<space|1fn><with|mode|text|in
      <with|mode|math|\<Omega\>>>,>>|<row|<cell|u>|<cell|=>|<cell|0<space|1em><with|mode|text|on
      <with|mode|math|\<partial\>\<Omega\>>>.>>>>
    </eqnarray*>
  </quotation>

  For <with|mode|math|n\<in\>\<bbb-N\>>, <with|mode|math|h=1/n> and
  <with|mode|math|\<cal-N\><rsub|h>=h\<bbb-Z\><rsup|2>\<cap\>\<Omega\>>. We
  set

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|1,h><rsup|+>u(x<rsub|1>,x<rsub|2>)>|<cell|\<assign\>>|<cell|<frac|1|h>(u(x<rsub|1>+h,x<rsub|2>)-u(x<rsub|1>,x<rsub|2>))>>|<row|<cell|\<partial\><rsub|1,h><rsup|->u(x<rsub|1>,x<rsub|2>)>|<cell|\<assign\>>|<cell|<frac|1|h>(u(x<rsub|1>-h,x<rsub|2>)-u(x<rsub|1>,x<rsub|2>))>>|<row|<cell|\<partial\><rsub|1,h><rsup|2>u(x<rsub|1>,x<rsub|2>)>|<cell|\<assign\>>|<cell|<frac|1|h<rsup|2>>(u(x<rsub|1>-h,x<rsub|2>)-2u(x<rsub|1>,x<rsub|2>)+u(x<rsub|1>+h,x<rsub|2>))>>|<row|<cell|-\<Delta\><rsub|h>u(x<rsub|1>,x<rsub|2>)>|<cell|\<assign\>>|<cell|<frac|1|h<rsup|2>>(4u(x<rsub|1>,x<rsub|2>)-u(x<rsub|1>-h,x<rsub|2>)-u(x<rsub|1>+h,x<rsub|2>)-u(x<rsub|1>,x<rsub|2>+h)-u(x<rsub|1>,x<rsub|2>-h))>>|<row|<cell|-\<Delta\><rsub|h>>|<cell|=>|<cell|<left|[><tabular*|<tformat|<table|<row|<cell|>|<cell|-1>|<cell|>>|<row|<cell|-1>|<cell|<phantom|->4>|<cell|-1>>|<row|<cell|>|<cell|-1>|<cell|>>>>><right|]><space|1fn><with|mode|text|(Finite
    difference star)>>>>>
  </eqnarray*>

  Find a grid function <with|mode|math|u:<wide|\<cal-N\><rsub|>|\<bar\>><rsub|h>\<assign\>\<bbb-Z\><rsup|2>\<cap\><wide|\<Omega\>|\<bar\>>\<rightarrow\>\<bbb-R\>>
  s.t. <with|mode|math|-\<Delta\><rsub|h>u<rsub|h>=f<rsub|h>> in
  <with|mode|math|\<cal-N\><rsub|h>>, <with|mode|math|u<rsub|h>=0> on
  <with|mode|math|\<partial\>\<cal-N\><rsub|h>=<wide|\<cal-N\>|\<bar\>><rsub|h>\<setminus\>\<cal-N\><rsub|h>>.
  By <em|lexicographic ordering> on <with|mode|math|\<cal-N\><rsub|h>={(i*h,j*h):i,j=1,\<ldots\>,n-1}>,
  <with|mode|math|u<rsub|h>\<in\>\<cal-N\><rsub|h>> is represented as
  <with|mode|math|\<b-u\>\<in\>\<bbb-R\><rsup|n>>
  (<with|mode|math|N=(n-1)<rsup|2>)>: <with|mode|math|\<b-u\>[i(n-1)+j]=u<rsub|h>(i*h,j*h)>.

  <\equation*>
    \<b-A\>\<b-u\>=\<b-f\> <with|mode|text|with>
    \<b-A\>=<matrix|<tformat|<table|<row|<cell|<phantom|->T>|<cell|-I>|<cell|>|<cell|>|<cell|>>|<row|<cell|-I>|<cell|<phantom|->T>|<cell|-I>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<ddots\>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|-I>|<cell|<phantom|->T>|<cell|-I>>|<row|<cell|>|<cell|>|<cell|>|<cell|-I>|<cell|<phantom|->T>>>>>\<in\>\<bbb-R\><rsup|N,N>\<backsimeq\>(\<bbb-R\><rsup|n-1,n-1>)<rsup|n-1,n-1>,
  </equation*>

  <\equation*>
    T=<matrix|<tformat|<table|<row|<cell|<phantom|->4>|<cell|-1>|<cell|>|<cell|>|<cell|>>|<row|<cell|-1>|<cell|<phantom|->4>|<cell|-1>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<ddots\>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|-1>|<cell|<phantom|->4>|<cell|-1>>|<row|<cell|>|<cell|>|<cell|>|<cell|-1>|<cell|<phantom|->4>>>>>\<in\>\<bbb-R\><rsup|n-1,n-1>,
  </equation*>

  <with|mode|math|I=diag(1,\<ldots\>,1)\<in\>\<bbb-R\><rsup|n-1,n-1>>.

  <\definition>
    <with|mode|math|\<b-A\>> sparse<with|mode|math|:\<Leftrightarrow\>><with|mode|math|\<b-A\>\<b-u\>>
    requires <with|mode|math|O(n)> operations.
  </definition>

  Here: max. 5 entries per row, bandwidth
  <with|mode|math|n><with|mode|math|\<Rightarrow\>><with|mode|math|\<b-A\>\<b-u\>>
  requires <with|mode|math|5n> operations. <with|mode|math|\<b-A\><rsup|-1>>
  requires <with|mode|math|O(N*n<rsup|2>)=O(N<rsup|2>)=O(n<rsup|4>)>
  operations and <with|mode|math|O(N*n)> storage (worse in 3D!).

  <subsection|Iterative Methods for Linear Systems>

  Let <with|mode|math|\<b-A\>\<in\>\<bbb-R\><rsup|N,N>> be regular,
  <with|mode|math|\<b-f\>\<in\>\<bbb-R\><rsup|N>> the right-hand side, and
  <with|mode|math|\<b-B\>\<in\>\<bbb-R\><rsup|N,N>> regular (the
  ``preconditioner''). For a given <with|mode|math|\<b-u\><rsup|0>\<in\>\<bbb-R\><rsup|N>>
  consider for <with|mode|math|k=1,2,\<ldots\>>
  <with|mode|math|\<b-u\><rsup|k>=\<b-u\><rsup|k-1>+\<b-B\>(\<b-f\>-\<b-A\>\<b-u\><rsup|k-1>)>.
  The error is, using <with|mode|math|\<b-u\>=\<b-A\><rsup|-1>\<b-f\>>,

  <\equation*>
    \<b-u\><rsup|k>-\<b-u\>=\<b-u\><rsup|k-1>-\<b-u\>+\<b-B\>(\<b-A\>\<b-u\>-\<b-A\>\<b-u\><rsup|k-1>)=(\<b-I\>-\<b-B\>\<b-A\>)(\<b-u\><rsup|k-1>-\<b-u\>)=(\<b-I\>-\<b-B\>\<b-A\>)<rsup|k>(\<b-u\><rsup|0>-\<b-u\>)
  </equation*>

  <\definition>
    <with|mode|math|\<sigma\>(\<b-A\>)\<assign\>{\<lambda\>\<in\>\<bbb-C\>:det(\<b-A\>-\<lambda\>\<b-I\>)=0}>,
    called the <em|spectrum>. <with|mode|math|\<rho\>(\<b-A\>)\<assign\>max{\|\<lambda\>\|:\<lambda\>\<in\>\<sigma\>(\<b-A\>)}>,
    called the spectral radius. <with|mode|math|\<kappa\>(\<b-A\>)\<assign\>\<rho\>(\<b-A\>)\<rho\>(\<b-A\><rsup|-1>)>,
    called the <em|condition number>. <with|mode|math|\<\|\|\>\<b-A\>\<\|\|\>\<assign\>sup<rsub|\<\|\|\>\<b-u\>\<\|\|\>=1>\<\|\|\>\<b-A\>\<b-u\>\<\|\|\>>,
    where <with|mode|math|\<\|\|\>\<b-u\>\<\|\|\>\<assign\><sqrt|\<b-u\><rsup|T>\<b-u\>>>,
    called the <em|spectral norm>.
  </definition>

  <\theorem>
    <with|mode|math|lim<rsub|k\<rightarrow\>\<infty\>><l2norm|(\<b-I\>-\<b-B\>\<b-A\>)<rsup|k>\<b-e\>||>=0>
    for any <with|mode|math|\<b-e\>\<in\>\<bbb-R\><rsup|N>><with|mode|math|\<Leftrightarrow\>><with|mode|math|\<rho\>(\<b-I\>-\<b-B\>\<b-A\>)\<less\>1>.
  </theorem>

  Application to <with|mode|math|\<b-B\>=\<theta\>\<b-I\>> with
  <with|mode|math|0\<less\>\<theta\>\<less\>1/\<rho\>(\<b-A\>)>:
  <with|mode|math|\<sigma\>(\<b-I\>-\<theta\>\<b-A\>)\<subset\>[-\<theta\>\<sigma\>(\<b-A\>),\<theta\>\<sigma\>(\<b-A\>)]><with|mode|math|\<Rightarrow\>>Iteration
  converges.

  <subsection|The Multigrid Idea>

  The basic assumption is that the fine-grid iteration has a hard time
  eliminating low-frequency errors. We consider the model problem. Define
  <with|mode|math|w<rsub|k><rsup|i>\<assign\><sqrt|h/2>*sin(h*k*i\<pi\>)>
  (<with|mode|math|i,k=1,\<ldots\>,n-1>, <with|mode|math|h=1/n>). Next,
  define <with|mode|math|\<b-w\><rsup|i,j>\<in\>\<bbb-R\><rsup|N>> by
  <with|mode|math|\<b-w\><rsup|i,j>[k(n-1)+l]=w<rsup|i><rsub|k>w<rsup|j><rsub|l>>.

  <\lemma>
    \;

    <\itemize>
      <item><with|mode|math|\<b-A\>\<b-w\><rsup|i,j>=\<lambda\><rsup|i,j>\<b-w\><rsup|i,j>>
      with <with|mode|math|\<lambda\><rsup|i,j>=4h<rsup|-2>(sin(i\<pi\>h/2)+sin<rsup|2>(j\<pi\>h/2))>,

      <item><with|mode|math|(\<b-w\><rsup|i,j>)> form a complete ONB,

      <item><with|mode|math|\<rho\>(\<b-A\>)=\<lambda\><rsub|max>=\<lambda\><rsup|n-1,n-1>=8h<rsup|-2>cos(\<pi\>h/2)=O(h<rsup|-2>)>,

      <with|mode|math|\<lambda\><rsub|min>=\<lambda\><rsup|1,1>=8h<rsup|-2>sin<rsup|2>(\<pi\>h/2)=O(1)>.
    </itemize>
  </lemma>

  Application to <with|mode|math|\<b-I\>-\<theta\>\<b-A\>>:
  <with|mode|math|\<sigma\>(\<b-I\>-\<theta\>\<b-A\>)\<subset\>[1-\<theta\>\<lambda\><rsub|max>,1-\<theta\>\<lambda\><rsub|min>]>,
  \ <with|mode|math|\<rho\>(\<b-I\>-\<theta\>\<b-A\>)=max{\|1-\<theta\>\<lambda\><rsub|max>\|,\|1-\<theta\>\<lambda\><rsub|min>\|}>.
  The optimal choice is

  <\equation*>
    \<theta\><rsup|\<ast\>>=<frac|2|\<lambda\><rsub|max>+\<lambda\><rsub|min>>=<frac|2|8h<rsup|-2>>=<frac|1|4>h<rsup|2>.
  </equation*>

  <with|mode|math|\<Rightarrow\>\<rho\>(\<b-I\>-\<theta\>\<b-A\>)=1-O(h<rsup|-2>)><with|mode|math|\<rightarrow\>>slow!

  <\lemma>
    Let <with|mode|math|\<cal-F\>={(i,j):max{i,j}\<geqslant\>n/2}>,
    <with|mode|math|X=span{\<b-w\><rsup|i,j>:(i,j)\<in\>\<cal-F\>}>, called
    the <em|subspace of strongly oscillating ``functions''>,
    <with|mode|math|dim X=3N/4>. Set <with|mode|math|\<theta\>=h<rsup|2>/8>,
    and assume for the initial error <with|mode|math|\<b-e\><rsup|0>=\<b-u\><rsup|0>-\<b-u\>\<in\>X>.
    Then, we have <with|mode|math|\<b-e\><rsup|k>=(\<b-I\>-\<theta\>\<b-A\>)\<b-e\><rsup|k-1>\<in\>X>
    and <with|mode|math|\<\|\|\>\<b-e\><rsup|k>\<\|\|\>\<leqslant\>3/4\<\|\|\>\<b-e\><rsup|k-1>\<\|\|\>\<leqslant\>(3/4)<rsup|k>\<\|\|\>\<b-e\><rsup|0>\<\|\|\>>.
  </lemma>

  <\proof>
    Consider

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<\|\|\>\<b-e\><rsup|k>\<\|\|\><rsup|2>>|<cell|=>|<cell|<big|sum><rsub|i,j>(\<b-e\><rsup|k>\<cdot\>\<b-w\><rsup|i,j>)<rsup|2>=<big|sum><rsub|i,j><left|[>(\<b-I\>-\<theta\>\<b-A\>)<rsup|k>\<b-e\><rsup|0>\<cdot\>\<b-w\><rsup|i,j><right|]><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|i,j><left|[>\<b-e\><rsup|0>\<cdot\>(\<b-I\>-\<theta\>\<b-A\>)<rsup|k>\<b-w\><rsup|i,j><right|]><rsup|2>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|(i,j)\<nin\>\<cal-F\>>(1-\<theta\>\<lambda\><rsup|i,j>)<rsup|2k><wide*|<left|[><wide*|\<b-e\><rsup|0>\<cdot\>\<b-w\><rsup|i,j>|\<wide-underbrace\>><rsub|\<b-e\><rsup|0>\<in\>X\<perp\>{\<b-w\><rsup|i,j>:i,j\<nin\>\<cal-F\>}
      ><right|]><rsup|2>|\<wide-underbrace\>><rsub|=0>+<big|sum><rsub|(i,j)\<in\>\<cal-F\>>(1-\<theta\>\<lambda\><rsup|i,j>)<rsup|2k><left|[>\<b-e\><rsup|0>\<cdot\>\<b-w\><rsup|i,j><right|]><rsup|2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<left|(>max<rsub|(i,j)\<in\>\<cal-F\>>(1-\<theta\>\<lambda\><rsup|i,j>)<rsup|2><right|)><rsup|k><l2norm|\<b-e\><rsup|0>||>.>>>>
    </eqnarray*>

    <\eqnarray*>
      <tformat|<table|<row|<cell|{\<lambda\><rsup|i,j>:(i,j)\<in\>\<cal-F\>}>|<cell|=>|<cell|4h<rsup|-2>{sin<rsup|2>(\<pi\>i*h/2)+sin<rsup|2>(\<pi\>j*h/2):i=1,\<ldots\>,n-1,j=n/2,\<ldots\>,n-1}>>|<row|<cell|>|<cell|\<subset\>>|<cell|4h<rsup|-2><left|[>sin<rsup|2>(\<pi\>/4),2]=h<rsup|-2>[2,8]>>>>
    </eqnarray*>

    <\eqnarray*>
      <tformat|<table|<row|<cell|max<rsub|(i,j)\<in\>\<cal-F\>>\|1-\<theta\>\<lambda\><rsup|i,j>\|>|<cell|=>|<cell|\|1-<frac|1|8>h<rsup|2>\<lambda\><rsup|1,n/2>\|>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|max<rsub|t\<in\>h<rsup|-2>[2,8]>\|1-<frac|1|8>h<rsup|2>t\|=<frac|3|4>.>>>>
    </eqnarray*>
  </proof>

  Consequence: Relaxation <with|mode|math|\<b-I\>-\<theta\>\<b-A\>>
  ``smoothing'' the error, slow for the smooth error components. Problem:
  <with|mode|math|X\<subset\>\<bbb-R\><rsup|N>> is not known in general.
  Idea: Continue relaxation with <em|coarse grid correction>.
</body>

<\initial>
  <\collection>
    <associate|language|american>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|1.1|?>>
    <associate|auto-3|<tuple|1.2|?>>
    <associate|auto-4|<tuple|1.3|?>>
    <associate|subsec:fd-model-problem|<tuple|1.1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|1.1<space|2spc>The Model Problem
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1.5fn>|1.2<space|2spc>Iterative Methods for
      Linear Systems <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1.5fn>|1.3<space|2spc>The Multigrid Idea
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>
    </associate>
  </collection>
</auxiliary>