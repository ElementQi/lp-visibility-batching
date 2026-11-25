import Mathlib.Analysis.Convex.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp

/-! Robust batched visibility in a shifted nonnegative convex image of a calibrated box. -/
namespace P1.BatchedVisibility

open scoped BigOperators

variable {n : ℕ}
abbrev Point (n : ℕ) := Fin n → ℝ

/-- The cap variables and their associated point in the shifted image. -/
def Feasible (K : Set (Point n)) (τ : ℝ) (u t : Point n) : Prop :=
  u ∈ K ∧ ∀ j, 0 ≤ t j ∧ t j ≤ τ ∧ t j ≤ u j

/-- Additive suboptimality, expressed without choosing a unique optimizer. -/
def Within (K : Set (Point n)) (τ : ℝ) (t : Point n) (gap : ℝ) : Prop :=
  ∀ u' t', Feasible K τ u' t' → ∑ j, t' j ≤ ∑ j, t j + gap

/-- A single feasible competitor controls the deficit at any witnessed head. -/
theorem cap_deficit {K : Set (Point n)} {τ gap : ℝ} {u t v : Point n} (i : Fin n)
    (conv : Convex ℝ K) (ht : Feasible K τ u t) (opt : Within K τ t gap)
    (hv : v ∈ K) (hv0 : ∀ j, 0 ≤ v j) (head : 1/4 ≤ v i) (small : τ ≤ 1) :
    (τ - t i) * (1/4 - (n : ℝ)*τ) ≤ gap := by
  classical
  let a := τ - t i
  have a0 : 0 ≤ a := sub_nonneg.mpr (ht.2 i).2.1
  have a1 : a ≤ 1 := by dsimp [a]; linarith [(ht.2 i).1]
  have b0 : 0 ≤ 1-a := sub_nonneg.mpr a1
  let z : Point n := (1-a) • u + a • v
  let q : Point n := fun j => (1-a)*t j + (if j=i then a/4 else 0)
  have zmem : z ∈ K := conv ht.1 hv b0 a0 (by ring)
  have feasible : Feasible K τ z q := by
    refine ⟨zmem, fun j => ?_⟩
    have base0 := mul_nonneg b0 (ht.2 j).1
    have base_le : (1-a)*t j ≤ t j := by nlinarith [mul_nonneg a0 (ht.2 j).1]
    have base_u := mul_le_mul_of_nonneg_left (ht.2 j).2.2 b0
    change 0 ≤ (1-a)*t j + _ ∧ (1-a)*t j + _ ≤ τ ∧
      (1-a)*t j + _ ≤ (1-a)*u j + a*v j
    by_cases h : j=i
    · subst j
      simp only [ite_true, eq_self]
      refine ⟨by linarith, ?_, ?_⟩
      · have ha : a = τ - t i := rfl
        linarith
      · nlinarith [mul_le_mul_of_nonneg_left head a0]
    · simp only [if_neg h, add_zero]
      exact ⟨base0, base_le.trans (ht.2 j).2.1, by linarith [mul_nonneg a0 (hv0 j)]⟩
  have sumq : ∑ j, q j = (1-a)*(∑ j, t j) + a/4 := by
    simp [q, Finset.sum_add_distrib, Finset.mul_sum]
  have sumt : ∑ j, t j ≤ (n : ℝ)*τ := by
    calc
      ∑ j, t j ≤ ∑ _j : Fin n, τ := Finset.sum_le_sum (fun j _ => (ht.2 j).2.1)
      _ = _ := by simp
  have improvement := opt z q feasible
  rw [sumq] at improvement
  change a * _ ≤ gap
  nlinarith [mul_le_mul_of_nonneg_left sumt a0]

/-- Exact optimal solutions saturate every witnessed slot simultaneously. -/
theorem exact_saturation {K : Set (Point n)} {τ : ℝ} {u t : Point n}
    (conv : Convex ℝ K) (ht : Feasible K τ u t) (opt : Within K τ t 0)
    (small : τ ≤ 1) (budget : (n : ℝ)*τ < 1/4) :
    ∀ i, (∃ v ∈ K, (∀ j, 0 ≤ v j) ∧ 1/4 ≤ v i) → t i = τ := by
  intro i ⟨v,hv,hv0,head⟩
  have h := cap_deficit i conv ht opt hv hv0 head small
  by_contra hne
  have deficit := lt_of_le_of_ne (ht.2 i).2.1 hne
  exact (not_lt_of_ge h) (mul_pos (sub_pos.mpr deficit) (sub_pos.mpr budget))

/-- A normalized opposite witness bounds the indicated nonnegative coordinate. -/
theorem opposite_bound (i : Fin n) {u p w : Point n} {d : ℝ}
    (hu : ∀ j, 0 ≤ u j) (hw : ∀ j, 0 ≤ w j ∧ w j ≤ 1)
    (head : 1/3 ≤ w i) (orth : ∑ j, w j*p j = 0)
    (err : ∀ j, u j - p j ≤ d) (hd : 0 ≤ d) : u i ≤ 3*(n : ℝ)*d := by
  have total : ∑ j, w j*u j ≤ (n : ℝ)*d := by
    calc
      ∑ j, w j*u j = ∑ j, (w j*p j + w j*(u j-p j)) := by
        apply Finset.sum_congr rfl
        intro j _
        ring
      _ = ∑ j, w j*(u j-p j) := by rw [Finset.sum_add_distrib, orth, zero_add]
      _ ≤ ∑ _j : Fin n, d := by
        apply Finset.sum_le_sum
        intro j _
        exact (mul_le_mul_of_nonneg_left (err j) (hw j).1).trans
          (by simpa using mul_le_mul_of_nonneg_right (hw j).2 hd)
      _ = _ := by simp
  have single : w i*u i ≤ ∑ j, w j*u j :=
    Finset.single_le_sum (f := fun j => w j*u j)
      (fun j _ => mul_nonneg (hw j).1 (hu j)) (Finset.mem_univ i)
  nlinarith [mul_le_mul_of_nonneg_right head (hu i)]

/-- The two promises of a price response, with the actual R182 cap and objective tolerance. -/
theorem price_labels {K : Set (Point n)} {τ d : ℝ} {u t : Point n}
    (hn : 0 < n) (cap : τ = 1/(32*(n : ℝ))) (conv : Convex ℝ K)
    (ht : Feasible K τ u t) (opt : Within K τ t (τ/32))
    (hd : 0 ≤ d) (precision : 3*(n : ℝ)*d < τ/2) (i : Fin n) :
    ((∃ v ∈ K, (∀ j, 0 ≤ v j) ∧ 1/4 ≤ v i) → τ/2 ≤ u i) ∧
    ((∃ p w : Point n, (∀ j, 0 ≤ w j ∧ w j ≤ 1) ∧ 1/3 ≤ w i ∧
      (∑ j, w j*p j) = 0 ∧ (∀ j, u j-p j ≤ d)) → u i < τ/2) := by
  have hn1 : (1:ℝ) ≤ n := by exact_mod_cast hn
  have hn0 : (n:ℝ) ≠ 0 := by linarith
  have budget : (n:ℝ)*τ = 1/32 := by rw [cap]; field_simp
  have small : τ ≤ 1 := by rw [cap]; apply (div_le_iff₀ (by positivity)).2; linarith
  have taupos : 0 < τ := by rw [cap]; positivity
  constructor
  · intro ⟨v,hv,hv0,head⟩
    have bound := cap_deficit i conv ht opt hv hv0 head small
    rw [budget] at bound
    linarith [(ht.2 i).2.2]
  · intro ⟨p,w,hw,head,orth,err⟩
    exact (opposite_bound i (fun j => (ht.2 j).1.trans (ht.2 j).2.2) hw head orth err hd).trans_lt precision

/-- On z=2x+y in the nonnegative 2-box, the unweighted sum uniquely favors x=0. -/
theorem uncapped_counterexample {x y z : ℝ} (hx : 0 ≤ x) (hy : y ≤ 2)
    (hz : z ≤ 2) (space : z=2*x+y) :
    (x+y+z ≤ 4) ∧ (4 ≤ x+y+z → x=0 ∧ y=2 ∧ z=2) := by
  constructor
  · linarith
  · intro h
    exact ⟨by linarith, by linarith, by linarith⟩


/-- The two slot blocks use a common coefficient vector. -/
def Boxed {m : ℕ} (l h : ℝ) (x : Point m) : Prop := ∀ j, l ≤ x j ∧ x j ≤ h

private theorem boxed_mix {m : ℕ} {x y : Point m} {l h a b : ℝ}
    (hx : Boxed l h x) (hy : Boxed l h y) (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a+b=1) :
    Boxed l h (a • x + b • y) := by
  intro j
  change l ≤ a*x j+b*y j ∧ a*x j+b*y j ≤ h
  constructor
  · calc
      l = a*l+b*l := by rw [← add_mul, hab, one_mul]
      _ ≤ _ := add_le_add (mul_le_mul_of_nonneg_left (hx j).1 ha)
        (mul_le_mul_of_nonneg_left (hy j).1 hb)
  · calc
      _ ≤ a*h+b*h := add_le_add (mul_le_mul_of_nonneg_left (hx j).2 ha)
        (mul_le_mul_of_nonneg_left (hy j).2 hb)
      _ = h := by rw [← add_mul, hab, one_mul]

def ShiftedBox {m : ℕ} (Bplus Bminus : Point m →ₗ[ℝ] Point n) (H ε : ℝ) :
    Set (Point n) := {u | ∃ c : Point m, Boxed (-H) H c ∧
      Boxed (-ε) 2 (Bplus c) ∧ Boxed (-ε) 2 (Bminus c) ∧ ∀ j, u j=Bplus c j+ε}

theorem shiftedBox_convex {m : ℕ} (Bplus Bminus : Point m →ₗ[ℝ] Point n) (H ε : ℝ) :
    Convex ℝ (ShiftedBox Bplus Bminus H ε) := by
  intro u hu v hv a b ha hb hab
  obtain ⟨c,hc,hpc,hmc,huc⟩ := hu
  obtain ⟨e,he,hpe,hme,hve⟩ := hv
  refine ⟨a • c+b • e, boxed_mix hc he ha hb hab, ?_, ?_, ?_⟩
  · simpa only [map_add, map_smul] using boxed_mix hpc hpe ha hb hab
  · simpa only [map_add, map_smul] using boxed_mix hmc hme ha hb hab
  · intro j
    simp only [map_add, map_smul, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
    rw [huc, hve]
    linarith only [show a*ε+b*ε=ε by rw [← add_mul, hab, one_mul]]

theorem shiftedBox_nonneg {m : ℕ} {Bplus Bminus : Point m →ₗ[ℝ] Point n} {H ε : ℝ}
    {u : Point n} (hu : u ∈ ShiftedBox Bplus Bminus H ε) : ∀ j, 0 ≤ u j := by
  obtain ⟨_,_,hp,_,shift⟩ := hu
  intro j
  rw [shift]; linarith [(hp j).1]

/-- Specialization to the actual paired coefficient box. The input facts correspond
to calibrated visibility, paired-space orthogonality and the entrywise error budget. -/
theorem calibrated_box_labels {m : ℕ} {Bplus Bminus : Point m →ₗ[ℝ] Point n}
    {H ε τ d : ℝ} {u t : Point n} (i : Fin n) (hn : 0 < n)
    (cap : τ=1/(32*(n:ℝ))) (eps : 0 ≤ ε)
    (ht : Feasible (ShiftedBox Bplus Bminus H ε) τ u t)
    (opt : Within (ShiftedBox Bplus Bminus H ε) τ t (τ/32))
    (hd : 0 ≤ d) (precision : 3*(n:ℝ)*d < τ/2) :
    ((∃ c, Boxed (-H) H c ∧ Boxed (-ε) 2 (Bplus c) ∧ Boxed (-ε) 2 (Bminus c) ∧
      1/4 ≤ Bplus c i) → τ/2 ≤ u i) ∧
    ((∃ p w : Point n, (∀ j, 0 ≤ w j ∧ w j ≤ 1) ∧ 1/3 ≤ w i ∧
      (∑ j, w j*p j)=0 ∧ (∀ j, u j-p j ≤ d)) → u i < τ/2) := by
  obtain ⟨pos,neg⟩ := price_labels hn cap (shiftedBox_convex Bplus Bminus H ε) ht opt hd precision i
  refine ⟨?_,neg⟩
  intro ⟨c,hc,hp,hm,head⟩
  let v : Point n := fun j => Bplus c j + ε
  have hv : v ∈ ShiftedBox Bplus Bminus H ε := ⟨c,hc,hp,hm,fun _ => rfl⟩
  exact pos ⟨v,hv,shiftedBox_nonneg hv,by dsimp [v]; linarith⟩


/-- Entrywise perturbations of a basis control the output on a bounded coefficient box. -/
theorem basis_error {m : ℕ} {A B : Fin m → Point n} {c : Point m} {H η : ℝ}
    (hc : ∀ k, |c k| ≤ H) (hη : 0 ≤ η) (close : ∀ k j, |B k j-A k j| ≤ η) (j : Fin n) :
    |(∑ k, c k*B k j) - (∑ k, c k*A k j)| ≤ (m:ℝ)*H*η := by
  rw [← Finset.sum_sub_distrib]
  calc
    _ ≤ ∑ k, |c k*B k j-c k*A k j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _k : Fin m, H*η := by
      apply Finset.sum_le_sum
      intro k _
      rw [← mul_sub, abs_mul]
      exact (mul_le_mul_of_nonneg_left (close k j) (abs_nonneg _)).trans
        (mul_le_mul_of_nonneg_right (hc k) hη)
    _ = _ := by simp [mul_assoc]

/-- Connect true-space positive and opposite witnesses to the batch LP in the approximate basis. -/
theorem labels_of_true_basis {m : ℕ} {Bplus Bminus Aplus Aminus : Point m →ₗ[ℝ] Point n}
    {H ε e τ : ℝ} {u t : Point n} (i : Fin n) (hn : 0 < n)
    (cap : τ=1/(32*(n:ℝ))) (he : 0 ≤ e) (heps : e ≤ ε) (tiny : e ≤ 1/12)
    (precision : 3*(n:ℝ)*(ε+e) < τ/2)
    (close : ∀ c, Boxed (-H) H c → ∀ j,
      |Bplus c j-Aplus c j| ≤ e ∧ |Bminus c j-Aminus c j| ≤ e)
    (ht : Feasible (ShiftedBox Bplus Bminus H ε) τ u t)
    (opt : Within (ShiftedBox Bplus Bminus H ε) τ t (τ/32)) :
    ((∃ c, Boxed (-H) H c ∧ Boxed 0 1 (Aplus c) ∧ Boxed 0 1 (Aminus c) ∧
      1/3 ≤ Aplus c i) → τ/2 ≤ u i) ∧
    ((∃ w : Point n, (∀ j, 0 ≤ w j ∧ w j ≤ 1) ∧ 1/3 ≤ w i ∧
      ∀ c, (∑ j, w j*Aplus c j)=0) → u i < τ/2) := by
  obtain ⟨pos,neg⟩ := calibrated_box_labels i hn cap (he.trans heps) ht opt
    (show 0 ≤ ε+e by linarith) precision
  constructor
  · intro ⟨c,hc,hp,hm,head⟩
    have box_transfer {x y : Point n} (hxy : ∀ j, |y j-x j| ≤ e) (hx : Boxed 0 1 x) :
        Boxed (-ε) 2 y := by
      intro j
      obtain ⟨lo,hi⟩ := abs_le.mp (hxy j)
      exact ⟨by linarith [(hx j).1],by linarith [(hx j).2]⟩
    apply pos
    refine ⟨c,hc,box_transfer (fun j => (close c hc j).1) hp,
      box_transfer (fun j => (close c hc j).2) hm,?_⟩
    linarith [(abs_le.mp (close c hc i).1).1]
  · intro ⟨w,hw,head,orth⟩
    obtain ⟨c,hc,_,_,shift⟩ := ht.1
    apply neg
    refine ⟨Aplus c,w,hw,head,orth c,fun j => ?_⟩
    rw [shift]
    linarith [(abs_le.mp (close c hc j).1).2]

/-- At the R182 cap budget every positively witnessed LP variable keeps six sevenths of its cap. -/
theorem robust_fraction {K : Set (Point n)} {τ : ℝ} {u t v : Point n} (i : Fin n)
    (conv : Convex ℝ K) (ht : Feasible K τ u t) (opt : Within K τ t (τ/32))
    (hv : v ∈ K) (hv0 : ∀ j, 0 ≤ v j) (head : 1/4 ≤ v i)
    (small : τ ≤ 1) (budget : (n:ℝ)*τ=1/32) : 6*τ/7 ≤ t i := by
  have h := cap_deficit i conv ht opt hv hv0 head small
  rw [budget] at h
  linarith


/-- The paper's concrete precision schedule supplies the scalar promises,
using its conservative allowance of at most 2n basis rows. -/
theorem manuscript_precision (n : ℕ) (hn : 0 < n) :
    let H : ℝ := 2^(100*n^2)
    let ε := 1/H^10
    let e := (2*(n:ℝ)*H)/H^100
    0 ≤ e ∧ e ≤ ε ∧ e ≤ 1/12 ∧ 3*(n:ℝ)*(ε+e) < (1/(32*(n:ℝ)))/2 := by
  dsimp only
  set H : ℝ := 2^(100*n^2)
  have hn1 : (1:ℝ) ≤ n := by exact_mod_cast hn
  have hpow : (n:ℝ) ≤ (2:ℝ)^n := by exact_mod_cast (Nat.lt_two_pow_self (n := n)).le
  have hexp : 10 + (n+n) ≤ 100*n^2 := by
    nlinarith [Nat.le_mul_of_pos_right n hn]
  have big : 1024*(n:ℝ)^2 ≤ H := by
    calc
      _ ≤ 1024*((2:ℝ)^n)^2 := by gcongr
      _ = (2:ℝ)^(10+(n+n)) := by rw [pow_add, pow_add, sq]; norm_num
      _ ≤ H := pow_le_pow_right₀ (by norm_num) hexp
  have nsq : (n:ℝ) ≤ (n:ℝ)^2 := by nlinarith [mul_nonneg (sub_nonneg.mpr hn1) (by positivity : 0 ≤ (n:ℝ))]
  have h1 : 1 ≤ H := by linarith
  have hp : 0 < H := by linarith
  have ep : 0 ≤ 2*(n:ℝ)*H / H^100 := by positivity
  have ee : 2*(n:ℝ)*H/H^100 ≤ 1/H^10 := by
    apply (div_le_div_iff₀ (by positivity) (by positivity)).2
    calc
      (2*(n:ℝ)*H)*H^10 ≤ (H*H)*H^10 := by gcongr; linarith
      _ = H^12 := by ring
      _ ≤ 1*H^100 := by rw [one_mul]; exact pow_le_pow_right₀ h1 (by decide)
  have inv : 1/H^10 ≤ 1/H := by
    apply one_div_le_one_div_of_le hp
    simpa using pow_le_pow_right₀ h1 (show 1 ≤ 10 by decide)
  refine ⟨ep,ee,?_,?_⟩
  · apply (ee.trans inv).trans
    apply (div_le_div_iff₀ hp (by norm_num)).2
    linarith
  · have upper : 3*(n:ℝ)*(1/H^10 + 2*(n:ℝ)*H/H^100) ≤ 6*n/H := by
      calc
        _ ≤ 3*(n:ℝ)*(1/H+1/H) := by exact mul_le_mul_of_nonneg_left (add_le_add inv (ee.trans inv)) (by positivity)
        _ = _ := by ring
    apply upper.trans_lt
    rw [div_div]
    apply (div_lt_div_iff₀ hp (by positivity)).2
    nlinarith


def rowEval {m : ℕ} (B : Fin m → Point n) : Point m →ₗ[ℝ] Point n where
  toFun c j := ∑ k, c k*B k j
  map_add' c d := by ext j; simp [add_mul, Finset.sum_add_distrib]
  map_smul' a c := by ext j; simp [Finset.mul_sum, mul_assoc]

/-- The batch decoder on explicit matrices with the published constants.
Both implications hold for one and the same feasible, approximately optimal LP answer. -/
theorem paper_box_decoder {m : ℕ} (i : Fin n) (hn : 0 < n) (hm : m ≤ 2*n)
    (Bp Bm Ap Am : Fin m → Point n) (u t : Point n) :
    let H : ℝ := 2^(100*n^2)
    let ε := 1/H^10
    let τ := 1/(32*(n:ℝ))
    (∀ k j, |Bp k j-Ap k j| ≤ 1/H^100 ∧ |Bm k j-Am k j| ≤ 1/H^100) →
    Feasible (ShiftedBox (rowEval Bp) (rowEval Bm) H ε) τ u t →
    Within (ShiftedBox (rowEval Bp) (rowEval Bm) H ε) τ t (τ/32) →
    ((∃ c, Boxed (-H) H c ∧ Boxed 0 1 (rowEval Ap c) ∧ Boxed 0 1 (rowEval Am c) ∧
      1/3 ≤ rowEval Ap c i) → τ/2 ≤ u i) ∧
    ((∃ w : Point n, (∀ j, 0 ≤ w j ∧ w j ≤ 1) ∧ 1/3 ≤ w i ∧
      ∀ c, (∑ j, w j*rowEval Ap c j)=0) → u i < τ/2) := by
  dsimp only
  set H : ℝ := 2^(100*n^2)
  intro close ht opt
  have precision := manuscript_precision n hn
  change 0 ≤ 2*(n:ℝ)*H/H^100 ∧ _ at precision
  obtain ⟨he,heps,htiny,hprec⟩ := precision
  apply labels_of_true_basis i hn rfl he heps htiny hprec _ ht opt
  intro c hc j
  have hc' : ∀ k, |c k| ≤ H := fun k => abs_le.mpr (hc k)
  have hp : 0 < H := by dsimp [H]; positivity
  have metaBound : (m:ℝ)*H*(1/H^100) ≤ 2*(n:ℝ)*H/H^100 := by
    rw [one_div, ← div_eq_mul_inv]
    apply div_le_div_of_nonneg_right _ (by positivity)
    apply mul_le_mul_of_nonneg_right _ hp.le
    exact_mod_cast hm
  constructor
  · exact (basis_error hc' (by positivity) (fun k l => (close k l).1) j).trans metaBound
  · exact (basis_error hc' (by positivity) (fun k l => (close k l).2) j).trans metaBound

end P1.BatchedVisibility

#print axioms P1.BatchedVisibility.price_labels
#print axioms P1.BatchedVisibility.exact_saturation

#print axioms P1.BatchedVisibility.calibrated_box_labels
#print axioms P1.BatchedVisibility.labels_of_true_basis
#print axioms P1.BatchedVisibility.robust_fraction
#print axioms P1.BatchedVisibility.manuscript_precision
#print axioms P1.BatchedVisibility.paper_box_decoder
