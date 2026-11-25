import Mathlib.Data.Real.Basic
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum

namespace P1.GraphCalibration

open scoped BigOperators

abbrev Slot (ι : Type*) := ι × Bool

def mirror {ι : Type*} (s : Slot ι) : Slot ι := (s.1, !s.2)

def signed {ι : Type*} (p : ι → ℝ) (s : Slot ι) : ℝ :=
  if s.2 then p s.1 else -p s.1

noncomputable def antisym {ι : Type*} (v : Slot ι → ℝ) (i : ι) : ℝ :=
  (v (i, true) - v (i, false)) / 2

theorem signed_antisym {ι : Type*} (v : Slot ι → ℝ) (s : Slot ι) :
    signed (antisym v) s = (v s - v (mirror s)) / 2 := by
  rcases s with ⟨i, h⟩
  cases h
  · simp [mirror, signed, antisym]
    ring
  · simp [mirror, signed, antisym]

theorem signed_mirror {ι : Type*} (p : ι → ℝ) (s : Slot ι) :
    signed p (mirror s) = -signed p s := by
  rcases s with ⟨i,h⟩
  cases h <;> simp [signed,mirror]

/-- A system of signed differences is feasible exactly when its doubled,
mirror-closed difference constraints have feasible potentials. -/
theorem signed_price_iff {ι E : Type*} (src dst : E → Slot ι) (cost : E → ℝ) :
    (∃ p : ι → ℝ, ∀ e, signed p (dst e) - signed p (src e) ≤ cost e) ↔
    (∃ v : Slot ι → ℝ, ∀ e, v (dst e) - v (src e) ≤ cost e ∧
      v (mirror (src e)) - v (mirror (dst e)) ≤ cost e) := by
  constructor
  · rintro ⟨p,h⟩
    refine ⟨signed p, fun e => ⟨h e, ?_⟩⟩
    rw [signed_mirror, signed_mirror]
    linarith [h e]
  · rintro ⟨v,h⟩
    refine ⟨antisym v, fun e => ?_⟩
    rw [signed_antisym, signed_antisym]
    linarith [(h e).1, (h e).2]

abbrev Form (d : ℕ) := Fin d → ℤ

def weight {d : ℕ} (a : Form d) : ℤ := ∑ i, |a i|

def eval {d : ℕ} (a : Form d) (x : Fin d → ℝ) : ℝ := ∑ i, (a i : ℝ) * x i

/-- A constant coordinate of value one can be included in x. Each cycle sum
of m integer affine edge forms is itself a test of weight at most m*A. -/
theorem cycle_type_transfer {d m : ℕ} {A : ℤ} (a : Fin m → Form d)
    (x y : Fin d → ℝ) (edges : ∀ j, weight (a j) ≤ A)
    (same : ∀ u : Form d, weight u ≤ (m:ℤ)*A → 0 ≤ eval u x → 0 ≤ eval u y)
    (cycle : 0 ≤ ∑ j, eval (a j) x) :
    0 ≤ ∑ j, eval (a j) y := by
  let u : Form d := fun i => ∑ j, a j i
  have hw : weight u ≤ (m:ℤ)*A := by
    calc
      weight u ≤ ∑ i, ∑ j, |a j i| :=
        Finset.sum_le_sum (fun i _ => Finset.abs_sum_le_sum_abs (fun j => a j i) _)
      _ = ∑ j, weight (a j) := Finset.sum_comm
      _ ≤ ∑ _j : Fin m, A := Finset.sum_le_sum (fun j _ => edges j)
      _ = (m:ℤ)*A := by simp
  have he (v : Fin d → ℝ) : eval u v = ∑ j, eval (a j) v := by
    simp only [eval,u,Int.cast_sum,Finset.sum_mul]
    exact Finset.sum_comm
  rw [← he]
  exact same u hw ((he x).symm ▸ cycle)

/-- A pivot selected by a proxy log gap larger than one increases physical
volume by at least 3/2 if the old-to-new normalized ratio has error ≤1/6. -/
theorem exchange_gain {old new proxy : ℝ} (hn : 0 < new)
    (gap : proxy ≤ 1/2) (close : old/new ≤ proxy+1/6) :
    3*old ≤ 2*new := by
  have h : old / new ≤ 2/3 := by linarith
  have := (div_le_iff₀ hn).mp h
  linarith


/-- Propagation of an entrywise error through one product in the Schur block.
The first approximate factor can have magnitude four and the second true
factor has magnitude at most three. -/
theorem product_error {x x' z z' h : ℝ} (hx : |x - x'| ≤ h)
    (hz : |z - z'| ≤ h) (a : |x'| ≤ 4) (b : |z| ≤ 3) :
    |x*z - x'*z'| ≤ 7*h := by
  calc
    |x*z - x'*z'| = |(x-x')*z + x'*(z-z')| := by ring_nf
    _ ≤ |x-x'| * |z| + |x'| * |z-z'| := by simpa only [abs_mul] using
      (abs_add_le ((x-x')*z) (x'*(z-z')))
    _ ≤ h*3 + 4*h :=
      add_le_add (mul_le_mul hx b (abs_nonneg _) (le_trans (abs_nonneg _) hx))
        (mul_le_mul a hz (abs_nonneg _) (by norm_num))
    _ = 7*h := by ring

theorem schur_error {k : ℕ} {y y' h : ℝ} (x x' z z' : Fin k → ℝ)
    (hy : |y-y'| ≤ h) (hx : ∀ i, |x i - x' i| ≤ h)
    (hz : ∀ i, |z i - z' i| ≤ h) (a : ∀ i, |x' i| ≤ 4) (b : ∀ i, |z i| ≤ 3) :
    |(y - ∑ i, x i*z i) - (y' - ∑ i, x' i*z' i)| ≤ (1+7*(k:ℝ))*h := by
  have hs : |∑ i, (x i*z i-x' i*z' i)| ≤ (k:ℝ)*(7*h) := by
    calc
      _ ≤ ∑ i, |x i*z i-x' i*z' i| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ _i : Fin k, (7*h) :=
        Finset.sum_le_sum (fun i _ => product_error (hx i) (hz i) (a i) (b i))
      _ = _ := by simp
  rw [Finset.sum_sub_distrib] at hs
  calc
    _ = |(y-y')-((∑ i, x i*z i)-(∑ i, x' i*z' i))| := by ring_nf
    _ ≤ |y-y'|+|((∑ i, x i*z i)-(∑ i, x' i*z' i))| := abs_sub _ _
    _ ≤ h + (k:ℝ)*(7*h) := add_le_add hy hs
    _ = _ := by ring

/-- Rowwise control of a perturbed identity solve, sufficient both for
injectivity (take b=0) and for controlling each recharted column. -/
theorem near_identity_solve {k : ℕ} (x : Fin k → ℝ) (f : Fin k → Fin k → ℝ)
    {b : ℝ} (hr : ∀ i, (∑ j, |f i j|) ≤ 1/2)
    (rhs : ∀ i, |x i + ∑ j, f i j*x j| ≤ b) (i : Fin k) : |x i| ≤ 2*b := by
  classical
  obtain ⟨m, _, hm⟩ := Finset.exists_max_image Finset.univ (fun j => |x j|)
    ⟨i, Finset.mem_univ _⟩
  have hsum : |∑ j, f m j*x j| ≤ 1/2 * |x m| := by
    calc
      _ ≤ ∑ j, |f m j*x j| := Finset.abs_sum_le_sum_abs _ _
      _ = ∑ j, |f m j| * |x j| := by simp only [abs_mul]
      _ ≤ ∑ j, |f m j| * |x m| := Finset.sum_le_sum (fun j hj =>
        mul_le_mul_of_nonneg_left (hm j hj) (abs_nonneg _))
      _ = (∑ j, |f m j|) * |x m| := by rw [Finset.sum_mul]
      _ ≤ _ := mul_le_mul_of_nonneg_right (hr m) (abs_nonneg _)
  have ht := abs_sub (x m + ∑ j, f m j*x j) (∑ j, f m j*x j)
  rw [add_sub_cancel_right] at ht
  have hi := hm i (Finset.mem_univ _)
  linarith [rhs m]

end P1.GraphCalibration

#print axioms P1.GraphCalibration.signed_price_iff
#print axioms P1.GraphCalibration.cycle_type_transfer
#print axioms P1.GraphCalibration.exchange_gain

#print axioms P1.GraphCalibration.schur_error
#print axioms P1.GraphCalibration.near_identity_solve
