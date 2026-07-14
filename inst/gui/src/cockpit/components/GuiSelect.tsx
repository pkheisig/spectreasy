import { Children, isValidElement, useEffect, useMemo, useRef, useState } from 'react'
import { Check, ChevronDown } from 'lucide-react'
import { createPortal } from 'react-dom'
import type { ChangeEvent, ReactNode, SelectHTMLAttributes } from 'react'

type Props = Omit<SelectHTMLAttributes<HTMLSelectElement>, 'children'> & { children: ReactNode }
type Option = { value: string; label: ReactNode; disabled: boolean }

export function GuiSelect({ children, className = '', value, defaultValue, disabled, onChange, 'aria-label': ariaLabel }: Props) {
  const buttonRef = useRef<HTMLButtonElement>(null)
  const [open, setOpen] = useState(false)
  const [internalValue, setInternalValue] = useState(String(defaultValue ?? ''))
  const [position, setPosition] = useState({ left: 0, top: 0, width: 180 })
  const options = useMemo<Option[]>(() => Children.toArray(children).flatMap((child) => {
    if (!isValidElement<{ value?: string | number; disabled?: boolean; children?: ReactNode }>(child) || child.type !== 'option') return []
    const optionValue = String(child.props.value ?? child.props.children ?? '')
    return [{ value: optionValue, label: child.props.children, disabled: Boolean(child.props.disabled) }]
  }), [children])
  const selectedValue = value == null ? internalValue : String(value)
  const selected = options.find((option) => option.value === selectedValue) ?? options[0]

  useEffect(() => {
    if (!open) return
    const close = (event: MouseEvent) => {
      if (!buttonRef.current?.contains(event.target as Node)) setOpen(false)
    }
    const escape = (event: KeyboardEvent) => {
      if (event.key === 'Escape') setOpen(false)
    }
    document.addEventListener('mousedown', close)
    document.addEventListener('keydown', escape)
    return () => {
      document.removeEventListener('mousedown', close)
      document.removeEventListener('keydown', escape)
    }
  }, [open])

  function toggle() {
    if (disabled) return
    const bounds = buttonRef.current?.getBoundingClientRect()
    if (bounds) setPosition({ left: bounds.left, top: bounds.bottom + 5, width: Math.max(bounds.width, 168) })
    setOpen((current) => !current)
  }

  function choose(option: Option) {
    if (option.disabled) return
    if (value == null) setInternalValue(option.value)
    onChange?.({ target: { value: option.value } } as ChangeEvent<HTMLSelectElement>)
    setOpen(false)
  }

  return (
    <>
      <button
        ref={buttonRef}
        type="button"
        className={`gui-select ${className} ${open ? 'is-open' : ''}`.trim()}
        aria-label={ariaLabel}
        aria-haspopup="listbox"
        aria-expanded={open}
        disabled={disabled}
        onClick={toggle}
      >
        <span>{selected?.label ?? selectedValue}</span>
        <ChevronDown size={15} />
      </button>
      {open && createPortal(
        <div className="gui-select-menu" role="listbox" style={position}>
          {options.map((option) => (
            <button
              type="button"
              role="option"
              aria-selected={option.value === selectedValue}
              className={option.value === selectedValue ? 'is-selected' : ''}
              disabled={option.disabled}
              key={option.value}
              onClick={() => choose(option)}
            >
              <span>{option.label}</span>
              {option.value === selectedValue && <Check size={15} />}
            </button>
          ))}
        </div>,
        document.body,
      )}
    </>
  )
}
