<TeXmacs|1.0.5>

<project|multigrid.tm>

<style|<tuple|article|skript|number-long-article>>

<\body>
  <assign|section-nr|5><section|Smoothing>

  Let <with|mode|math|V<rsub|j>=V=<big|sum><rsub|n=1><rsup|N>W<rsub|n>>, e.g.
  <with|mode|math|W<rsub|n>=span \<varphi\><rsub|n>> with norm
  <with|mode|math|<l2norm|\<cdot\>||>>. Further, consider

  <\eqnarray*>
    <tformat|<table|<row|<cell|a:V\<times\>V>|<cell|\<rightarrow\>>|<cell|\<bbb-R\>,<with|mode|text|
    symmetric positive definite>, with <enorm|v||>=<sqrt|a(v,v)>.>>|<row|<cell|A:V>|<cell|\<rightarrow\>>|<cell|V,<with|mode|text|
    defined by ><ip|A<rsub|n>v|w||>=a(v,w)>>|<row|<cell|A<rsub|n>:W<rsub|n>>|<cell|\<rightarrow\>>|<cell|W<rsub|n>,<with|mode|text|
    defined by ><ip|A<rsub|n>v<rsub|n>|w<rsub|n>||>=a(v<rsub|n>,w<rsub|n>)>>|<row|<cell|Q<rsub|n>:V>|<cell|\<rightarrow\>>|<cell|W<rsub|n>,<with|mode|text|
    defined by ><ip|Q*v|w<rsub|n>||>=<ip|v|w<rsub|n>||>>>|<row|<cell|P<rsub|n>:V>|<cell|\<rightarrow\>>|<cell|W<rsub|n>,
    <with|mode|text| defined by ><ip|P<rsub|n>v|w<rsub|n>||>=a(v,w<rsub|n>)>>>>
  </eqnarray*>

  <subsection|Method of succesive subspace correction>

  For <with|mode|math|u<rprime|'>\<in\>V> and
  <with|mode|math|l\<in\>V<rprime|'>> compute for
  <with|mode|math|n=1,\<ldots\>,N> the vector
  <with|mode|math|c<rsub|n>\<in\>W<rsub|n>>, such that

  <\eqnarray*>
    <tformat|<table|<row|<cell|a(c<rsub|n>,w<rsub|n>)>|<cell|=>|<cell|l(w<rsub|n>)-a(u<rsup|n-1>,w<rsub|n>)>>|<row|<cell|<with|mode|text|and
    then >u<rsup|n>>|<cell|\<assign\>>|<cell|u<rsup|n-1>+c<rsub|n>>>>>
  </eqnarray*>

  In operator form, we compute <with|mode|math|r<rsup|0>\<in\>V> with
  <with|mode|math|<ip|r<rsup|0>|v||>=l(v)-a(u<rsup|0>,v)>. Then, for
  <with|mode|math|k=1,2,3,\<ldots\>>, compute

  <\eqnarray*>
    <tformat|<table|<row|<cell|c<rsub|n>>|<cell|\<assign\>>|<cell|A<rsup|-1><rsub|n>Q<rsub|n>r<rsup|n-1>>>|<row|<cell|c<rsup|n>>|<cell|\<assign\>>|<cell|c<rsup|n+1>+c<rsub|n>>>|<row|<cell|r<rsup|n>>|<cell|\<assign\>>|<cell|r<rsup|n-1>-A*c<rsub|n>>>>>
  </eqnarray*>

  Let <with|mode|math|u\<in\>V> be the solution of
  <with|mode|math|a(u,v)=l(v)> for all <with|mode|math|v\<in\>V>. We obtain

  <\equation*>
    u<rsup|k>-u=(id-R*A)(u<rsup|k-1>-u)
  </equation*>

  with

  <\equation*>
    (id-R*A)=(id-P<rsub|N>)\<cdots\>(id-P<rsub|1>),
  </equation*>

  with <with|mode|math|P<rsub|i>\<assign\>A<rsub|n><rsup|-1>Q<rsub|n>A>.

  <\theorem>
    <label|the:smoothing-prop-gen-smoother>Assume

    <\enumerate-alpha>
      <item><with|mode|math|\<forall\>w\<in\>V<space|0.6spc>\<exists\>w<rsub|n>\<in\>W<rsub|n>:<space|0.6spc>w=<big|sum><rsub|n=1><rsup|N>w<rsub|n>>
      and <with|mode|math|<big|sum><rsub|n=1><rsup|N><l2norm|w<rsub|n>||2>\<leqslant\>C<rsub|0><l2norm|w||2>=C<rsub|0><l2norm|<big|sum><rsub|n>W<rsub|n>||2>>,

      <item><with|mode|math|<l2norm|\<frak-X\>|2|>\<leqslant\>C<rsub|1>>,
      where <with|mode|math|\<frak-X\>\<in\>\<bbb-R\><rsup|N\<times\>N>> is
      the so-called <em|interaction matrix>, with

      <\equation*>
        \<frak-X\><rsub|n,m>=<choice|<tformat|<table|<row|<cell|0>|<cell|P<rsub|n>P<rsub|m>=0,>>|<row|<cell|1>|<cell|<with|mode|text|otherwise>.>>>>>
      </equation*>
    </enumerate-alpha>

    Then, the smoothing property

    <\equation*>
      <enorm|(id-R*A)||2>\<leqslant\>a((id-\<omega\>\<lambda\><rsup|-1>A)v,v)<space|1em><with|mode|text|with><space|1em>\<omega\>=<frac|1|C<rsub|0>C<rsub|1><rsup|2>>\<in\>(0,1)<space|1em><with|mode|text|and><space|1em>\<lambda\>=<l2norm|A||>
    </equation*>

    holds.
  </theorem>

  <\proof>
    <em|Step 1.> Sei <with|mode|math|E<rsub|0>\<assign\>id>,
    <with|mode|math|E<rsub|1>\<assign\>id-P<rsub|1>>,
    <with|mode|math|E<rsub|m>=(1-P<rsub|m>)E<rsub|n><with|color|red|???>\<cdots\>>
    and <with|mode|math|id-R*A=E<rsub|N>>. Then,

    <\equation*>
      P<rsub|n>E<rsub|m>=0,<space|1em>E<rsub|m-1>-E<rsub|m>=(id-(id-P<rsub|n>))E<rsub|m-1>=P<rsub|m>E<rsub|m-1>
    </equation*>

    and

    <\equation>
      id-E<rsub|n>=<big|sum><rsub|m=1><rsup|n>E<rsub|m-1>-E<rsub|m>=<big|sum><rsub|m=1><rsup|n>P<rsub|m>E<rsub|m-1><label|eq:smoothing-id-e>
    </equation>

    Now,

    <\equation*>
      a(P<rsub|n>E<rsub|n-1>v,E<rsub|n-1>v)=a(E<rsub|m-1>v,E<rsub|n-1>v)-a(E<rsub|m>v,E<rsub|n-1>v)
    </equation*>

    and

    <\eqnarray*>
      <tformat|<table|<row|<cell|a(E<rsub|n>v,E<rsub|n-1>v)>|<cell|=>|<cell|a(E<rsub|n>v,E<rsub|n-1>v)-<wide*|a(<wide*|P<rsub|n>E<rsub|m>|\<wide-underbrace\>><rsub|=0>,E<rsub|n-1>v)|\<wide-underbrace\>><rsub|=a(E<rsub|n>v,P<rsub|n>E<rsub|n-1>v)>>>|<row|<cell|>|<cell|=>|<cell|a(E<rsub|n>v,E<rsub|n-1>v-P<rsub|m>E<rsub|n-1>v)=a(E<rsub|n>v,E<rsub|n>v)=a(E<rsub|n>v,P<rsub|n>E<rsub|n-1>v).>>>>
    </eqnarray*>

    As a consequence,

    <\eqnarray*>
      <tformat|<table|<row|<cell|a(v,v)-a(E<rsub|N>v,E<rsub|N>v)>|<cell|=>|<cell|<big|sum><rsub|n-1><rsup|N>a(E<rsub|n-1>v,E<rsub|n-1>v)-a(E<rsub|n>v,E<rsub|n>v)>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|n=1><rsup|N>a(P<rsub|n>E<rsub|n-1>v,E<rsub|n-1>v).>>>>
    </eqnarray*>

    <em|Step 2>. Show <with|mode|math|<l2norm|A*v||2>\<leqslant\>C<rsub|0>C<rsub|1><rsup|2>\<lambda\><big|sum><rsub|n=1><rsup|N>a(P<rsub|n>E<rsub|n-1>v,E<rsub|n-1>v)>.

    For <with|mode|math|v\<in\>V>, we set <with|mode|math|w=A*v> and choose
    <with|mode|math|w<rsub|n>\<in\>W<rsub|n>> with
    <with|mode|math|<big|sum><rsub|n=1><rsup|N><l2norm|w<rsub|n>||2>\<leqslant\>C<rsub|0><l2norm|w||2>>,
    where <with|mode|math|w=<big|sum><rsub|n=1><rsup|N>w<rsub|n>>.

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|w||2>>|<cell|=>|<cell|<l2norm|A*v||2>=<ip|A*v|A*v||>=a(v,w)=<big|sum><rsub|n=1><rsup|N>a(v,w<rsub|n>)>>|<row|<cell|>|<cell|<below|=|(<reference|eq:smoothing-id-e>)>>|<cell|<big|sum><rsub|n=1><rsup|N>a<left|(><left|(>E<rsub|n+1>+<big|sum><rsub|m=1><rsup|n-1>P<rsub|m>E<rsub|m-1><right|)>v,w<rsub|n><right|)>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|n=1><rsup|N><left|(><wide*|a(E<rsub|n-1>v,w<rsub|n>)|\<wide-underbrace\>><rsub|=a(P<rsub|n>E<rsub|n-1>v,w<rsub|n>)>+<big|sum><rsub|m=1><rsup|n-1>a(P<rsub|m>E<rsub|m-1>v,w<rsub|n>)<right|)>>>|<row|<cell|>|<cell|=>|<cell|<big|sum><rsub|n=1><rsup|N><big|sum><rsub|m=1><rsup|n>a(P<rsub|m>E<rsub|m-1>v,w<rsub|n>)\<frak-X\><rsub|n,m>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<big|sum><rsub|n=1><rsup|N><big|sum><rsub|m=1><rsup|n><wide*|<enorm|P<rsub|m>E<rsub|m-1>v||>|\<wide-underbrace\>><rsub|x<rsub|m>>*<wide*|<enorm|w<rsub|n>||>|\<wide-underbrace\>><rsub|y<rsub|m>>\<frak-X\><rsub|n,m>
      >>|<row|<cell|>|<cell|<below|\<leqslant\>|(\<ast\>)>>|<cell|<l2norm|\<frak-X\>|2|><sqrt|<big|sum><rsub|n=1><rsup|N><wide*|<enorm|P<rsub|m>E<rsub|m-1>||2>|\<wide-underbrace\>><rsub|=a(P<rsub|m>,E<rsub|m-1>v)>><sqrt|<big|sum><wide*|<enorm|w<rsub|n>||2>|\<wide-underbrace\>><rsub|=<ip|A*w<rsub|n>|w<rsub|n>||>\<leqslant\>\<lambda\><l2norm|w<rsub|n>||2>>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|C<rsub|1><sqrt|<big|sum>a(P<rsub|n>E<rsub|n-1>v,E<rsub|n>v><sqrt|\<lambda\>C<rsub|0><l2norm|w||2>>,>>>>
    </eqnarray*>

    where the step <with|mode|math|(\<ast\>)> uses, for
    <with|mode|math|x,y\<in\>\<bbb-R\><rsup|n>>,

    <\equation*>
      <big|sum><big|sum>x<rsub|n>\<frak-X\><rsub|n,m>y<rsub|m>=x<rsup|T>\<frak-X\>y\<leqslant\><l2norm|x|2|><l2norm|\<frak-X\>y|2|>\<leqslant\><l2norm|x|2|><l2norm|y|2|><l2norm|\<frak-X\>|2|>.
    </equation*>

    Altogether, step 2 is proven.

    <em|Step 3.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|<enorm|v||2>-<enorm|E<rsub|N><with|color|red|??>||2>>|<cell|<below|=|<with|mode|text|Step
      1>>>|<cell|<big|sum>a(P<rsub|n>E<rsub|n-1>v,E<rsub|n-1>v)\<geqslant\><frac|1|C<rsub|0>C<rsub|1><rsup|2>\<lambda\>><l2norm|A*v||2>=<frac|\<omega\>|\<lambda\>>.>>>>
    </eqnarray*>

    So,

    <\equation*>
      <enorm|id-R*A||2>=<enorm|E<rsub|N>v||2>\<leqslant\><enorm|v||2>-<frac|\<omega\>|\<lambda\>><l2norm|A*v||2>=a<left|(><left|(>id-<frac|\<omega\>|\<lambda\>>A<right|)>v,v<right|)>.
    </equation*>

    \;
  </proof>

  <subsection|Parallel Subspace Correction>

  This chapter discusses the so-called <em|BPX method>. (BPX for Bramble, P?,
  Xu?) For a given residual <with|mode|math|r<rsup|k>>, compute in parallel
  for <with|mode|math|n=1,\<ldots\>,N> the correction
  <with|mode|math|c<rsub|n>=A<rsup|-1>Q<rsub|n>r<rsup|k>>. Then, set

  <\equation*>
    u<rsup|k+1>=u<rsup|k>+c<space|1em><with|mode|text|with><space|1em>c=\<gamma\><big|sum>c<rsub|n>,
  </equation*>

  where <with|mode|math|\<gamma\>\<in\>(0,1)> is a damping factor. The new
  residual is <with|mode|math|r<rsup|k+1>=r<rsup|k>-A*c>. So,

  <\equation*>
    (u<rsup|k+1>-u)=(id-<wide*|R|\<wide-underbrace\>><rsub|=R<rsub|par>>A)(u<rsup|k>-u)
  </equation*>

  with

  <\equation*>
    R=\<gamma\><big|sum>A<rsub|n><rsup|-1>Q<rsub|n>A*A<rsup|-1>=\<gamma\><big|sum>P<rsub|m>A<rsup|-1>.
  </equation*>

  For comparison,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<with|mode|text|Successive Subspace
    Correction:><space|1em>id-R*A>|<cell|=>|<cell|<big|prod>(id-P<rsub|n>),>>|<row|<cell|<with|mode|text|Parallel
    Subspace Correction:><space|1em>id-R*A>|<cell|=>|<cell|id-\<gamma\><big|sum>P<rsub|n>.>>>>
  </eqnarray*>

  However, it turns out that the parallel correction only achieves
  convergence rates which are roughly half of those achieved by successive
  correction.

  <\theorem>
    Assume a), b) of Theorem <reference|the:smoothing-prop-gen-smoother>.
    Choose <with|mode|math|\<gamma\>> such that
    <with|mode|math|C<rsub|1>\<gamma\>\<in\>(0,2)>. Then, the smoothing
    property of SSC holds for PSC with <with|mode|math|\<omega\>=<frac|\<gamma\>|C<rsub|0>>(2-\<gamma\>C<rsub|1>)>.
  </theorem>

  <\proof>
    <em|Step 1:> Show\ 

    <\equation*>
      <l2norm|v||2>\<leqslant\><wide*|C<rsub|0>|\<wide-underbrace\>><rsub|<with|mode|text|cf.
      a)>><big|sum><rsub|n=1><rsup|N><l2norm|Q<rsub|n>v||2>.
    </equation*>

    Choose <with|mode|math|v=<big|sum>w<rsub|n>>. With a),

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|v||2>>|<cell|=>|<cell|<ip|v|<big|sum>w<rsub|n>||>=<big|sum><ip|v|w<rsub|n>||>=<big|sum><ip|v|Q<rsub|n>w<rsub|n>||>=<big|sum><ip|Q<rsub|n>v|w<rsub|n>||>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<sqrt|<big|sum><l2norm|Q<rsub|n>v||2>><wide*|<sqrt|<big|sum><l2norm|w<rsub|n>||2>>|\<wide-underbrace\>><rsub|\<leqslant\>C<rsub|0><l2norm|w||2>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<sqrt|<big|sum><l2norm|Q<rsub|n>v||2>><sqrt|C<rsub|0><l2norm|w||2>>.>>>>
    </eqnarray*>

    So,

    <\equation*>
      <l2norm|v||2>\<leqslant\>C<rsub|0><big|sum><l2norm|Q<rsub|n>w<rsub|n>||>*<l2norm|v||>.<value|huh>
    </equation*>

    <em|Step 2.> We deduce

    <\eqnarray*>
      <tformat|<table|<row|<cell|<l2norm|v||2>>|<cell|\<leqslant\>>|<cell|C<rsub|0><big|sum><l2norm|Q<rsub|n>v||2>=C<rsub|0><big|sum><l2norm|A<rsub|n><rsup|1/2>A<rsub|n><rsup|-1/2>Q<rsub|n>v||2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|C<rsub|0>\<lambda\><big|sum><l2norm|A<rsub|n><rsup|-1/2>Q<rsub|n>v||2>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|C<rsub|0>\<lambda\><big|sum><ip|A<rsup|-1><rsub|n>Q<rsub|n>v|Q<rsub|n>v||>\<leqslant\>C<rsub|0><frac|\<lambda\>|\<gamma\>><ip|R*v|v||>.>>>>
    </eqnarray*>

    using the fact that

    <\equation*>
      <l2norm|A<rsub|n><rsup|1/2>||2>=sup<rsub|w<rsub|n>\<in\>W<rsub|n>,<l2norm|w<rsub|n>||>=1><ip|A<rsub|n>w<rsub|n>|w<rsub|n>||>\<leqslant\>sup<rsub|<l2norm|v||>=1><ip|A*v|v||>=\<lambda\>.
    </equation*>

    <em|Step 3.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|a(R*v,R*v)>|<cell|=>|<cell|\<gamma\><rsup|2><big|sum><rsub|n=1><rsup|N><big|sum><rsub|m=1><rsup|N>a(<wide*|P<rsub|n>P<rsub|m>|\<wide-underbrace\>><rsub|<with|mode|text|could
      be added because of <with|mode|math|\<frak-X\>>>>A<rsup|-1><rsub|n>Q<rsub|n>v,A<rsup|-1><rsub|m>Q<rsub|m>v)\<frak-X\><rsub|n,m>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|\<gamma\><rsup|2><big|sum><rsub|n=1><rsup|N><big|sum><rsub|m=1><rsup|N><enorm|A<rsup|-1><rsub|m>Q<rsub|m>v||>*<enorm|A<rsup|-1><rsub|n>Q<rsub|n>v||>\<frak-X\><rsub|n,m>>>|<row|<cell|>|<cell|<below|\<leqslant\>|<with|mode|text|CSU>>>|<cell|\<gamma\><rsup|2><l2norm|\<frak-X\>|2|><sqrt|<big|sum><rsub|m><enorm|A<rsub|m><rsup|-1>Q<rsub|m>v||2>><sqrt|<big|sum><rsub|n><enorm|A<rsup|-1><rsub|n>Q<rsub|n>v||2>>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|\<gamma\><rsup|2>C<rsub|1><wide*|a(A<rsup|-1><rsub|m>Q<rsub|m>v,A<rsup|-1><rsub|n>Q<rsub|n>v)|\<wide-underbrace\>><rsub|<ip|A<rsup|-1><rsub|n>Q<rsub|n>v|v||>>>>|<row|<cell|>|<cell|=>|<cell|\<gamma\>C<rsub|1><ip|R*v|v||>.>>>>
    </eqnarray*>

    <em|Step 4.>

    <\eqnarray*>
      <tformat|<table|<row|<cell|2<ip|R*v|v||>-a(R*v,R*v)>|<cell|\<geqslant\>>|<cell|(2-C<rsub|1>\<gamma\>)<ip|R*v|v||>>>|<row|<cell|>|<cell|\<geqslant\>>|<cell|(2-C<rsub|1>\<gamma\>)<frac|\<gamma\>|C<rsub|0>\<lambda\>><l2norm|v||2>>>|<row|<cell|>|<cell|=>|<cell|<frac|\<omega\>|\<lambda\>><l2norm|v||2>.>>>>
    </eqnarray*>

    So,

    <\eqnarray*>
      <tformat|<table|<row|<cell|<enorm|(id-R*A*)v||2>>|<cell|<below|=|v\<rightarrow\>A*v>>|<cell|<enorm|v||2>-2<ip|R*A*v|A*v||>+a(R*A*v,R*A*v)>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<enorm|v||2>-<frac|\<omega\>|\<lambda\>><l2norm|A*v||2>>>|<row|<cell|>|<cell|=>|<cell|a<left|(><left|(>id-<frac|\<omega\>|\<lambda\>>A<right|)>v,v<right|)>.>>>>
    </eqnarray*>

    (I.e. we're better than Richardson.(?))
  </proof>

  <\lemma>
    On a uniform mesh in <with|mode|math|\<bbb-R\><rsup|2>>, we have for a
    space <with|mode|math|V<rsub|h>> of linear finite elements

    <\eqnarray*>
      <tformat|<table|<row|<cell|C<rsub|2><big|sum><rsub|z\<in\>N<rsub|C>>h<rsup|2>v(z)<rsup|2>\<leqslant\>>|<cell|<l2norm|v|0|2>>|<cell|\<leqslant\>C<rsub|3><big|sum>h<rsup|2>v(z)<rsup|2>>>>>
    </eqnarray*>

    for all <with|mode|math|v\<in\>V<rsub|h>>.
  </lemma>

  <\proof>
    On an element <with|mode|math|\<Omega\><rsub|C>>, we have

    <\equation*>
      <l2norm|v|\<Omega\><rsub|C>|2>\<approx\>\|\<Omega\><rsub|C>\|<big|sum><rsub|z\<in\>N<rsub|C>>v(z)<rsup|2>\<approx\>h<rsup|2><big|sum><rsub|z\<in\>N<rsub|C>>v(z)<rsup|2>.
    </equation*>

    \;
  </proof>

  <\corollary>
    For the Gauÿ-Seidel method with <with|mode|math|w<rsub|n>=span
    \<varphi\><rsub|n>>, we have for <with|mode|math|w<rsub|n>=v(z<rsub|n>)\<varphi\><rsub|n>>
    (and thus <with|mode|math|v=<big|sum>v(z<rsub|n>)\<varphi\><rsub|n>>)

    <\equation*>
      <big|sum><l2norm|w<rsub|n>||2>\<leqslant\>C<rsub|0><l2norm|v||2>
    </equation*>
  </corollary>

  <\proof>
    We have <with|mode|math|<l2norm|w<rsub|n>||2>=v(z)<rsup|2><l2norm|\<varphi\><rsub|n>||2>\<leqslant\>C<rsub|4>h<rsup|2>v(z)<rsup|2>>uniformly,
    depending on connectivity. So, <with|mode|math|C<rsub|0>=C<rsub|2>C<rsub|4>>.
  </proof>

  \;
</body>

<\initial>
  <\collection>
    <associate|language|american>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|6|?>>
    <associate|auto-2|<tuple|6.1|?>>
    <associate|auto-3|<tuple|6.2|?>>
    <associate|eq:smoothing-id-e|<tuple|6.1|?>>
    <associate|the:smoothing-prop-gen-smoother|<tuple|6.1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Smoothing>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|6.1<space|2spc>Method of succesive
      subspace correction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1.5fn>|6.2<space|2spc>Parallel Subspace
      Correction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>
    </associate>
  </collection>
</auxiliary>