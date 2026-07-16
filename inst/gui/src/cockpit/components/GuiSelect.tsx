import { Children, isValidElement, useEffect, useId, useMemo, useRef, useState } from 'react'
import { Check, ChevronDown } from 'lucide-react'
import { createPortal } from 'react-dom'
import type { ChangeEvent, KeyboardEvent as ReactKeyboardEvent, ReactNode, SelectHTMLAttributes } from 'react'

type Props = Omit<SelectHTMLAttributes<HTMLSelectElement>, 'children'> & { children: ReactNode }
type Option = { value: string; label: ReactNode; disabled: boolean }

export function GuiSelect({ children, className = '', value, defaultValue, disabled, onChange, 'aria-label': ariaLabel }: Props) {
  const buttonRef = useRef<HTMLButtonElement>(null)
  const menuRef = useRef<HTMLDivElement>(null)
  const optionRefs = useRef<Array<HTMLButtonElement | null>>([])
  const menuId = useId()
  const [open, setOpen] = useState(false)
  const [activeIndex, setActiveIndex] = useState(0)
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
      const target = event.target as Node
      if (!buttonRef.current?.contains(target) && !menuRef.current?.contains(target)) setOpen(false)
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

  useEffect(() => {
    if (open) optionRefs.current[activeIndex]?.focus()
  }, [activeIndex, open])

  function positionMenu() {
    const bounds = buttonRef.current?.getBoundingClientRect()
    if (bounds) setPosition({ left: bounds.left, top: bounds.bottom + 5, width: Math.max(bounds.width, 168) })
  }

  function openMenu(preferredIndex = options.findIndex((option) => option.value === selectedValue)) {
    if (disabled || options.length === 0) return
    positionMenu()
    setActiveIndex(Math.max(0, preferredIndex))
    setOpen(true)
  }

  function toggle() {
    if (disabled) return
    if (open) setOpen(false)
    else openMenu()
  }

  function choose(option: Option) {
    if (option.disabled) return
    if (value == null) setInternalValue(option.value)
    const target = { value: option.value } as HTMLSelectElement
    onChange?.({ target, currentTarget: target } as ChangeEvent<HTMLSelectElement>)
    setOpen(false)
    window.requestAnimationFrame(() => buttonRef.current?.focus())
  }

  function moveActive(delta: number) {
    if (!options.length) return
    let next = activeIndex
    do next = (next + delta + options.length) % options.length
    while (options[next]?.disabled && next !== activeIndex)
    setActiveIndex(next)
  }

  function handleButtonKeyDown(event: ReactKeyboardEvent<HTMLButtonElement>) {
    if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
      event.preventDefault()
      const selectedIndex = Math.max(0, options.findIndex((option) => option.value === selectedValue))
      openMenu(event.key === 'ArrowDown' ? selectedIndex : Math.max(0, options.length - 1))
    }
  }

  function handleMenuKeyDown(event: ReactKeyboardEvent<HTMLDivElement>) {
    if (event.key === 'Escape') {
      event.preventDefault()
      setOpen(false)
      buttonRef.current?.focus()
    } else if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
      event.preventDefault()
      moveActive(event.key === 'ArrowDown' ? 1 : -1)
    } else if (event.key === 'Home' || event.key === 'End') {
      event.preventDefault()
      setActiveIndex(event.key === 'Home' ? 0 : options.length - 1)
    } else if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault()
      const option = options[activeIndex]
      if (option) choose(option)
    }
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
        aria-controls={open ? menuId : undefined}
        disabled={disabled}
        onClick={toggle}
        onKeyDown={handleButtonKeyDown}
      >
        <span>{selected?.label ?? selectedValue}</span>
        <ChevronDown size={15} />
      </button>
      {open && createPortal(
        <div id={menuId} ref={menuRef} className="gui-select-menu" role="listbox" style={position} onKeyDown={handleMenuKeyDown}>
          {options.map((option, index) => (
            <button
              ref={(element) => { optionRefs.current[index] = element }}
              type="button"
              role="option"
              aria-selected={option.value === selectedValue}
              className={option.value === selectedValue ? 'is-selected' : ''}
              disabled={option.disabled}
              tabIndex={index === activeIndex ? 0 : -1}
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
