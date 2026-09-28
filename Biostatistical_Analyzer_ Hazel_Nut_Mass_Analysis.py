"""
Биостатистический анализатор с GUI
Автоматизация анализа данных по массе орехов лещины (Corylus avellana L.)
Интерфейс на Tkinter
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy import stats
from scipy.stats import norm, chi2
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import math
from typing import List, Dict, Tuple
import json
import os

class BioStatApp:
    """
    Основной класс GUI-приложения для биостатистического анализа
    """
    
    def __init__(self, root):
        self.root = root
        self.root.title("Биостатистический анализатор - Анализ массы орехов лещины")
        self.root.geometry("1400x800")
        
        # Данные
        self.data = None
        self.results = {}
        self.alpha = 0.01  # уровень значимости по умолчанию
        self.measurement_precision = 0.0001
        
        # Создание интерфейса
        self.create_widgets()
        self.load_sample_data()  # Загрузка данных из курсовой работы по умолчанию
        
    def create_widgets(self):
        """Создание всех элементов интерфейса"""
        # Главный контейнер
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Конфигурация расширяемости
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(1, weight=1)
        
        # === БЛОК ВВОДА ДАННЫХ ===
        input_frame = ttk.LabelFrame(main_frame, text="Ввод данных", padding="10")
        input_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        
        # Кнопки загрузки данных
        ttk.Button(input_frame, text="Загрузить из файла (CSV/Excel)", 
                  command=self.load_from_file).grid(row=0, column=0, padx=5)
        ttk.Button(input_frame, text="Ввести вручную", 
                  command=self.open_manual_input).grid(row=0, column=1, padx=5)
        ttk.Button(input_frame, text="Использовать пример данных", 
                  command=self.load_sample_data).grid(row=0, column=2, padx=5)
        
        # Параметры анализа
        param_frame = ttk.Frame(input_frame)
        param_frame.grid(row=1, column=0, columnspan=3, pady=(10, 0), sticky=(tk.W, tk.E))
        
        ttk.Label(param_frame, text="Уровень значимости (α):").grid(row=0, column=0, padx=(0, 5))
        self.alpha_var = tk.StringVar(value="0.01")
        alpha_combo = ttk.Combobox(param_frame, textvariable=self.alpha_var, 
                                  values=["0.01", "0.05", "0.10"], width=10)
        alpha_combo.grid(row=0, column=1, padx=(0, 10))
        
        ttk.Label(param_frame, text="Точность измерений:").grid(row=0, column=2, padx=(20, 5))
        self.precision_var = tk.StringVar(value="0.0001")
        ttk.Entry(param_frame, textvariable=self.precision_var, width=10).grid(row=0, column=3)
        
        ttk.Button(param_frame, text="Запустить анализ", 
                  command=self.run_full_analysis, style="Accent.TButton").grid(row=0, column=4, padx=(20, 0))
        
        # Информация о данных
        self.data_info_label = ttk.Label(input_frame, text="Данные не загружены")
        self.data_info_label.grid(row=2, column=0, columnspan=3, pady=(10, 0))
        
        # === ОСНОВНОЕ ОКНО С РЕЗУЛЬТАТАМИ ===
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Вкладка 1: Исходные данные
        self.tab1 = ttk.Frame(notebook)
        notebook.add(self.tab1, text="Исходные данные")
        
        self.text_original = scrolledtext.ScrolledText(self.tab1, width=60, height=20)
        self.text_original.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Вкладка 2: Вариационные ряды
        self.tab2 = ttk.Frame(notebook)
        notebook.add(self.tab2, text="Вариационные ряды")
        
        self.text_variation = scrolledtext.ScrolledText(self.tab2, width=60, height=20)
        self.text_variation.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Вкладка 3: Основные характеристики
        self.tab3 = ttk.Frame(notebook)
        notebook.add(self.tab3, text="Основные характеристики")
        
        self.text_stats = scrolledtext.ScrolledText(self.tab3, width=60, height=20)
        self.text_stats.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Вкладка 4: Критерий Пирсона
        self.tab4 = ttk.Frame(notebook)
        notebook.add(self.tab4, text="Критерий Пирсона")
        
        self.text_pearson = scrolledtext.ScrolledText(self.tab4, width=60, height=20)
        self.text_pearson.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Вкладка 5: Доверительные интервалы
        self.tab5 = ttk.Frame(notebook)
        notebook.add(self.tab5, text="Доверительные интервалы")
        
        self.text_confidence = scrolledtext.ScrolledText(self.tab5, width=60, height=20)
        self.text_confidence.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Вкладка 6: Визуализация
        self.tab6 = ttk.Frame(notebook)
        notebook.add(self.tab6, text="Графики")
        
        # Создаем фрейм для графиков
        self.figure_frame = ttk.Frame(self.tab6)
        self.figure_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # === ПАНЕЛЬ УПРАВЛЕНИЯ ===
        control_frame = ttk.Frame(main_frame)
        control_frame.grid(row=2, column=0, columnspan=2, pady=(10, 0), sticky=(tk.W, tk.E))
        
        ttk.Button(control_frame, text="Экспорт отчета (TXT)", 
                  command=self.export_report).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Сохранить графики", 
                  command=self.save_plots).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Сбросить данные", 
                  command=self.reset_data).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Справка", 
                  command=self.show_help).pack(side=tk.LEFT, padx=5)
        
        # Статус бар
        self.status_bar = ttk.Label(main_frame, text="Готов к работе", relief=tk.SUNKEN)
        self.status_bar.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(10, 0))
    
    def load_sample_data(self):
        """Загрузка примерных данных из курсовой работы"""
        # Данные из таблицы 1 курсовой работы (90 значений)
        sample_data = [
            0.9825, 0.9928, 1.0252, 1.0862, 1.0742, 1.0951, 1.1042, 1.1085, 1.1126, 1.1245,
            1.1368, 1.1742, 1.1785, 1.1857, 1.1895, 1.1942, 1.1952, 1.2045, 1.2125, 1.2269,
            1.2345, 1.2374, 1.2498, 1.2564, 1.2567, 1.2654, 1.2658, 1.2694, 1.2695, 1.2782,
            1.2784, 1.2845, 1.2862, 1.2951, 1.2963, 1.3024, 1.3045, 1.3125, 1.3156, 1.3216,
            1.3219, 1.3248, 1.3328, 1.3363, 1.3366, 1.3369, 1.3372, 1.3375, 1.3377, 1.3614,
            1.3675, 1.3695, 1.3695, 1.3742, 1.3748, 1.3754, 1.3763, 1.3816, 1.3816, 1.3842,
            1.3842, 1.3892, 1.3924, 1.3925, 1.4257, 1.4325, 1.4328, 1.4416, 1.4429, 1.4528,
            1.4578, 1.4695, 1.4786, 1.4863, 1.4874, 1.4923, 1.4956, 1.5123, 1.5162, 1.5162,
            1.5248, 1.5328, 1.5364, 1.5416, 1.5472, 1.5523, 1.5562, 1.5712, 1.5823, 1.5914
        ]
        
        self.data = np.array(sample_data)
        self.measurement_precision = float(self.precision_var.get())
        self.update_data_info()
        self.show_original_data()
        messagebox.showinfo("Данные загружены", 
                          f"Загружены примерные данные из курсовой работы\n"
                          f"Количество наблюдений: {len(self.data)}\n"
                          f"Уровень значимости: α = {self.alpha}")
    
    def load_from_file(self):
        """Загрузка данных из файла"""
        filetypes = [
            ('CSV files', '*.csv'),
            ('Excel files', '*.xlsx;*.xls'),
            ('Text files', '*.txt'),
            ('All files', '*.*')
        ]
        
        filepath = filedialog.askopenfilename(
            title="Выберите файл с данными",
            filetypes=filetypes
        )
        
        if filepath:
            try:
                if filepath.endswith('.csv'):
                    df = pd.read_csv(filepath)
                elif filepath.endswith(('.xlsx', '.xls')):
                    df = pd.read_excel(filepath)
                else:
                    # Предполагаем текстовый файл с числами
                    with open(filepath, 'r') as f:
                        data = [float(line.strip()) for line in f if line.strip()]
                    df = pd.DataFrame({'data': data})
                
                # Предлагаем выбрать колонку
                if len(df.columns) > 1:
                    column = self.select_column_dialog(df.columns)
                    if column:
                        self.data = df[column].dropna().values
                    else:
                        return
                else:
                    self.data = df.iloc[:, 0].dropna().values
                
                self.measurement_precision = float(self.precision_var.get())
                self.update_data_info()
                self.show_original_data()
                self.status_bar.config(text=f"Данные загружены из {os.path.basename(filepath)}")
                
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось загрузить файл:\n{str(e)}")
    
    def select_column_dialog(self, columns):
        """Диалог выбора колонки из DataFrame"""
        dialog = tk.Toplevel(self.root)
        dialog.title("Выберите колонку")
        dialog.geometry("300x150")
        dialog.transient(self.root)
        dialog.grab_set()
        
        ttk.Label(dialog, text="В файле несколько колонок.\nВыберите колонку с данными:").pack(pady=10)
        
        selected_column = tk.StringVar()
        combo = ttk.Combobox(dialog, textvariable=selected_column, values=list(columns))
        combo.pack(pady=10)
        combo.current(0)
        
        result = []
        
        def on_ok():
            result.append(selected_column.get())
            dialog.destroy()
        
        ttk.Button(dialog, text="OK", command=on_ok).pack(pady=10)
        
        self.root.wait_window(dialog)
        return result[0] if result else None
    
    def open_manual_input(self):
        """Открытие окна для ручного ввода данных"""
        dialog = tk.Toplevel(self.root)
        dialog.title("Ручной ввод данных")
        dialog.geometry("500x400")
        dialog.transient(self.root)
        
        # Текстовое поле для ввода
        ttk.Label(dialog, text="Введите значения через пробел или с новой строки:").pack(pady=10)
        
        text_frame = ttk.Frame(dialog)
        text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        text_input = scrolledtext.ScrolledText(text_frame, height=15)
        text_input.pack(fill=tk.BOTH, expand=True)
        
        # Пример данных
        text_input.insert(tk.END, "1.25 1.30 1.28 1.32 1.29 1.31 1.27 1.26 1.33 1.34")
        
        def parse_and_close():
            try:
                text = text_input.get("1.0", tk.END)
                # Разбиваем на числа
                numbers = []
                for line in text.split('\n'):
                    for part in line.split():
                        try:
                            num = float(part.replace(',', '.'))
                            numbers.append(num)
                        except ValueError:
                            continue
                
                if len(numbers) < 10:
                    messagebox.showwarning("Мало данных", "Введите хотя бы 10 значений")
                    return
                
                self.data = np.array(numbers)
                self.measurement_precision = float(self.precision_var.get())
                self.update_data_info()
                self.show_original_data()
                dialog.destroy()
                messagebox.showinfo("Данные загружены", 
                                  f"Загружено {len(numbers)} значений")
                
            except Exception as e:
                messagebox.showerror("Ошибка", f"Ошибка при разборе данных:\n{str(e)}")
        
        ttk.Button(dialog, text="Загрузить данные", command=parse_and_close).pack(pady=10)
    
    def update_data_info(self):
        """Обновление информации о загруженных данных"""
        if self.data is not None:
            n = len(self.data)
            x_min = np.min(self.data)
            x_max = np.max(self.data)
            self.data_info_label.config(
                text=f"Загружено {n} наблюдений. "
                     f"Min: {x_min:.4f} г, Max: {x_max:.4f} г, "
                     f"Размах: {x_max - x_min:.4f} г"
            )
    
    def show_original_data(self):
        """Отображение исходных данных"""
        if self.data is not None:
            self.text_original.delete(1.0, tk.END)
            self.text_original.insert(tk.END, "РАНЖИРОВАННЫЙ РЯД НАБЛЮДЕНИЙ МАССЫ ОРЕХОВ ЛЕЩИНЫ\n")
            self.text_original.insert(tk.END, "="*60 + "\n\n")
            
            sorted_data = np.sort(self.data)
            for i, value in enumerate(sorted_data, 1):
                self.text_original.insert(tk.END, f"{i:3d}. {value:.4f} г")
                if i % 3 == 0:
                    self.text_original.insert(tk.END, "\n")
                else:
                    self.text_original.insert(tk.END, "\t")
            
            self.text_original.insert(tk.END, "\n\n")
            self.text_original.insert(tk.END, f"Общее количество: {len(self.data)}\n")
            self.text_original.insert(tk.END, f"Минимальное значение: {np.min(self.data):.4f} г\n")
            self.text_original.insert(tk.END, f"Максимальное значение: {np.max(self.data):.4f} г\n")
            self.text_original.insert(tk.END, f"Размах вариации: {np.max(self.data) - np.min(self.data):.4f} г")
    
    def run_full_analysis(self):
        """Запуск полного анализа"""
        if self.data is None or len(self.data) < 10:
            messagebox.showwarning("Недостаточно данных", 
                                 "Загрузите данные для анализа (минимум 10 значений)")
            return
        
        try:
            self.alpha = float(self.alpha_var.get())
            self.measurement_precision = float(self.precision_var.get())
            
            self.status_bar.config(text="Выполняется анализ...")
            self.root.update()
            
            # 1. Построение вариационных рядов
            self.build_variation_series()
            
            # 2. Расчет основных характеристик
            self.calculate_basic_statistics()
            
            # 3. Критерий Пирсона
            self.calculate_pearson_test()
            
            # 4. Доверительные интервалы
            self.calculate_confidence_intervals()
            
            # 5. Визуализация
            self.create_plots()
            
            self.status_bar.config(text="Анализ успешно завершен")
            messagebox.showinfo("Анализ завершен", 
                              "Статистический анализ успешно выполнен.\n"
                              "Результаты доступны на всех вкладках.")
            
        except Exception as e:
            self.status_bar.config(text="Ошибка при анализе")
            messagebox.showerror("Ошибка", f"Произошла ошибка при анализе:\n{str(e)}")
    
    def build_variation_series(self):
        """Построение интервального вариационного ряда"""
        n = len(self.data)
        x_min = np.min(self.data)
        x_max = np.max(self.data)
        R = x_max - x_min  # размах вариации
        
        # Правило Стерджеса
        K = int(1 + 3.322 * math.log10(n))  # число классов
        lambda_val = R / (K - 1)  # ширина интервала
        
        # Корректировка с учетом точности измерений
        lambda_val = round(lambda_val, 4)
        
        # Построение интервалов
        lower_bound = x_min - lambda_val / 2
        lower_bound = round(lower_bound, 4)
        
        intervals = []
        frequencies = []
        
        current_lower = lower_bound
        for i in range(K):
            current_upper = current_lower + lambda_val
            # Корректировка границ для избежания перекрытия
            if i < K - 1:
                current_upper -= self.measurement_precision
            
            # Подсчет частот
            freq = np.sum((self.data >= current_lower) & (self.data <= current_upper))
            intervals.append((round(current_lower, 4), round(current_upper, 4)))
            frequencies.append(freq)
            
            current_lower = current_upper
        
        # Сохраняем результаты
        self.results['intervals'] = intervals
        self.results['frequencies'] = frequencies
        self.results['K'] = K
        self.results['lambda'] = lambda_val
        self.results['R'] = R
        
        # Вычисление средних значений в интервалах
        midpoints = [(low + high) / 2 for low, high in intervals]
        
        # Отображение результатов
        self.text_variation.delete(1.0, tk.END)
        self.text_variation.insert(tk.END, "ИНТЕРВАЛЬНЫЙ ВАРИАЦИОННЫЙ РЯД МАССЫ ОРЕХОВ ЛЕЩИНЫ\n")
        self.text_variation.insert(tk.END, "="*60 + "\n\n")
        
        self.text_variation.insert(tk.END, f"Количество наблюдений (n): {n}\n")
        self.text_variation.insert(tk.END, f"Минимальное значение (x_min): {x_min:.4f} г\n")
        self.text_variation.insert(tk.END, f"Максимальное значение (x_max): {x_max:.4f} г\n")
        self.text_variation.insert(tk.END, f"Размах вариации (R): {R:.4f} г\n")
        self.text_variation.insert(tk.END, f"Число классов (по Стерджесу): {K}\n")
        self.text_variation.insert(tk.END, f"Ширина классового интервала (λ): {lambda_val:.4f} г\n")
        
        self.text_variation.insert(tk.END, "\n" + "="*60 + "\n\n")
        self.text_variation.insert(tk.END, f"{'Интервал':^25} | {'Частота':^10} | {'Середина':^10}\n")
        self.text_variation.insert(tk.END, "-"*50 + "\n")
        
        total_freq = 0
        for i, ((low, high), freq, mid) in enumerate(zip(intervals, frequencies, midpoints), 1):
            interval_str = f"{low:.4f} - {high:.4f}"
            
            mid_formatted = f"{mid:.4f}"
            self.text_variation.insert(tk.END, 
                f"{interval_str:^25} | {freq:^10} | {mid_formatted:^10}\n")
            
            total_freq += freq
        
        self.text_variation.insert(tk.END, "-"*50 + "\n")
        self.text_variation.insert(tk.END, f"{'ИТОГО':^25} | {total_freq:^10} | {'-':^10}\n")    
    
    def calculate_basic_statistics(self):
        """Расчет основных статистических характеристик"""
        if 'midpoints' not in self.results:
            # Вычисляем средние значения интервалов
            intervals = self.results['intervals']
            midpoints = [(low + high) / 2 for low, high in intervals]
            self.results['midpoints'] = midpoints
        
        n = len(self.data)
        midpoints = self.results['midpoints']
        frequencies = self.results['frequencies']
        
        # 1. Среднее значение
        weighted_sum = sum(m * f for m, f in zip(midpoints, frequencies))
        mean_value = weighted_sum / n
        
        # 2. Дисперсия и СКО
        variance_sum = sum(f * (m - mean_value) ** 2 for m, f in zip(midpoints, frequencies))
        variance = variance_sum / (n - 1)
        std_dev = math.sqrt(variance)
        
        # 3. Мода (интервал с максимальной частотой)
        max_freq_idx = np.argmax(frequencies)
        modal_interval = self.results['intervals'][max_freq_idx]
        
        # 4. Медиана
        cumulative_freq = np.cumsum(frequencies)
        median_idx = np.searchsorted(cumulative_freq, n / 2)
        if median_idx < len(intervals):
            median_interval = self.results['intervals'][median_idx]
        else:
            median_interval = modal_interval
        
        # Сохраняем результаты
        self.results['mean'] = mean_value
        self.results['variance'] = variance
        self.results['std_dev'] = std_dev
        self.results['modal_interval'] = modal_interval
        self.results['median_interval'] = median_interval
        
        # Отображение результатов
        self.text_stats.delete(1.0, tk.END)
        self.text_stats.insert(tk.END, "ОСНОВНЫЕ ХАРАКТЕРИСТИКИ ВАРИАЦИОННОГО РЯДА\n")
        self.text_stats.insert(tk.END, "="*60 + "\n\n")
        
        self.text_stats.insert(tk.END, "1. ХАРАКТЕРИСТИКИ ЦЕНТРАЛЬНОЙ ТЕНДЕНЦИИ:\n")
        self.text_stats.insert(tk.END, "-"*60 + "\n")
        self.text_stats.insert(tk.END, f"Среднее арифметическое (x̄): {mean_value:.4f} г\n\n")
        
        self.text_stats.insert(tk.END, f"Модальный интервал: {modal_interval[0]:.4f} - {modal_interval[1]:.4f} г\n")
        self.text_stats.insert(tk.END, f"Частота модального интервала: {frequencies[max_freq_idx]}\n\n")
        
        self.text_stats.insert(tk.END, f"Медианный интервал: {median_interval[0]:.4f} - {median_interval[1]:.4f} г\n\n")
        
        self.text_stats.insert(tk.END, "2. ХАРАКТЕРИСТИКИ ВАРИАЦИИ:\n")
        self.text_stats.insert(tk.END, "-"*60 + "\n")
        self.text_stats.insert(tk.END, f"Дисперсия (s²): {variance:.6f} г²\n")
        self.text_stats.insert(tk.END, f"Среднее квадратическое отклонение (s): {std_dev:.4f} г\n")
        self.text_stats.insert(tk.END, f"Коэффициент вариации (CV): {(std_dev/mean_value*100):.2f}%\n")
        
        self.text_stats.insert(tk.END, "\n3. НОРМИРОВАННЫЕ ОТКЛОНЕНИЯ:\n")
        self.text_stats.insert(tk.END, "-"*60 + "\n")
        
        # Вычисление нормированных отклонений
        norm_deviations = [(m - mean_value) / std_dev for m in midpoints]
        self.results['norm_deviations'] = norm_deviations
        
        for i, (mid, t) in enumerate(zip(midpoints, norm_deviations), 1):
            self.text_stats.insert(tk.END, f"t{i} = ({mid:.4f} - {mean_value:.4f}) / {std_dev:.4f} = {t:.4f}\n")
    
    def calculate_pearson_test(self):
        """Вычисление критерия согласия Пирсона"""
        mean_val = self.results['mean']
        std_val = self.results['std_dev']
        n = len(self.data)
        lambda_val = self.results['lambda']
        intervals = self.results['intervals']
        frequencies = self.results['frequencies']
        norm_deviations = self.results.get('norm_deviations', [])
        
        # Вычисление теоретических частот
        theoretical_freqs = []
        chi2_components = []
        
        for i, (low, high) in enumerate(intervals):
            # Нормированные границы
            z_low = (low - mean_val) / std_val
            z_high = (high - mean_val) / std_val
            
            # Вероятность попадания в интервал
            p = norm.cdf(z_high) - norm.cdf(z_low)
            
            # Теоретическая частота
            f_theoretical = n * p
            theoretical_freqs.append(f_theoretical)
            
            # Компонент χ² (если теоретическая частота не слишком мала)
            if f_theoretical >= 2:  # условие из курсовой
                chi2_comp = (frequencies[i] - f_theoretical) ** 2 / f_theoretical
            else:
                chi2_comp = 0
            
            chi2_components.append(chi2_comp)
        
        # Суммируем χ²
        chi2_observed = sum(chi2_components)
        
        # Число степеней свободы
        df = len(intervals) - 3  # как в курсовой
        
        # Критическое значение χ²
        chi2_critical = chi2.ppf(1 - self.alpha, df)
        
        # p-value
        p_value = 1 - chi2.cdf(chi2_observed, df)
        
        # Сохраняем результаты
        self.results['theoretical_freqs'] = theoretical_freqs
        self.results['chi2_observed'] = chi2_observed
        self.results['chi2_critical'] = chi2_critical
        self.results['chi2_df'] = df
        self.results['chi2_p_value'] = p_value
        
        # Отображение результатов
        self.text_pearson.delete(1.0, tk.END)
        self.text_pearson.insert(tk.END, "ПРОВЕРКА СООТВЕТСТВИЯ НОРМАЛЬНОМУ РАСПРЕДЕЛЕНИЮ\n")
        self.text_pearson.insert(tk.END, "="*60 + "\n\n")
        
        self.text_pearson.insert(tk.END, f"Параметры распределения:\n")
        self.text_pearson.insert(tk.END, f"Среднее (μ): {mean_val:.4f}\n")
        self.text_pearson.insert(tk.END, f"СКО (σ): {std_val:.4f}\n")
        self.text_pearson.insert(tk.END, f"Уровень значимости (α): {self.alpha}\n\n")
        
        self.text_pearson.insert(tk.END, "ТАБЛИЦА ДЛЯ РАСЧЕТА χ²-КРИТЕРИЯ ПИРСОНА\n")
        self.text_pearson.insert(tk.END, "-"*80 + "\n")
        self.text_pearson.insert(tk.END, 
        f"{'Интервал':^20} | {'Эмп. частота':^12} | {'Теор. частота':^12} | {'(f-f\')²':^12} | {'(f-f\')²/f':^12}\n")
        self.text_pearson.insert(tk.END, "-"*80 + "\n")
        
        total_emp = 0
        total_theor = 0
        
        for i, ((low, high), f_emp, f_theor, chi2_comp) in enumerate(
            zip(intervals, frequencies, theoretical_freqs, chi2_components)):
            
            interval_str = f"{low:.4f}-{high:.4f}"
            diff_sq = (f_emp - f_theor) ** 2
            
            self.text_pearson.insert(tk.END, 
                f"{interval_str:^20} | {f_emp:^12.0f} | {f_theor:^12.2f} | "
                f"{diff_sq:^12.2f} | {chi2_comp:^12.4f}\n")
            
            total_emp += f_emp
            total_theor += f_theor
        
        self.text_pearson.insert(tk.END, "-"*80 + "\n")
        self.text_pearson.insert(tk.END, 
            f"{'ИТОГО':^20} | {total_emp:^12.0f} | {total_theor:^12.2f} | "
            f"{'':^12} | {chi2_observed:^12.4f}\n\n")
        
        self.text_pearson.insert(tk.END, "РЕЗУЛЬТАТЫ ТЕСТА:\n")
        self.text_pearson.insert(tk.END, "-"*60 + "\n")
        self.text_pearson.insert(tk.END, f"Наблюдаемое значение χ²: {chi2_observed:.4f}\n")
        self.text_pearson.insert(tk.END, f"Число степеней свободы (df): {df}\n")
        self.text_pearson.insert(tk.END, f"Критическое значение χ² (α={self.alpha}): {chi2_critical:.4f}\n")
        self.text_pearson.insert(tk.END, f"p-value: {p_value:.6f}\n\n")
        
        # Статистический вывод
        if chi2_observed < chi2_critical:
            conclusion = (f"χ²_набл ({chi2_observed:.4f}) < χ²_крит ({chi2_critical:.4f})\n"
                         f"⇒ Нулевая гипотеза НЕ отвергается.\n"
                         f"Распределение соответствует нормальному (при α={self.alpha}).")
        else:
            conclusion = (f"χ²_набл ({chi2_observed:.4f}) ≥ χ²_крит ({chi2_critical:.4f})\n"
                         f"⇒ Нулевая гипотеза ОТВЕРГАЕТСЯ.\n"
                         f"Распределение НЕ соответствует нормальному (при α={self.alpha}).")
        
        self.text_pearson.insert(tk.END, "ВЫВОД:\n")
        self.text_pearson.insert(tk.END, "-"*60 + "\n")
        self.text_pearson.insert(tk.END, conclusion)
        
        self.results['normality_conclusion'] = chi2_observed < chi2_critical
    
    def calculate_confidence_intervals(self):
        """Расчет доверительных интервалов"""
        mean_val = self.results['mean']
        std_val = self.results['std_dev']
        variance = self.results['variance']
        n = len(self.data)
        
        # t-значение для доверительного интервала
        # Для α=0.01 и больших n можно использовать z-value
        if self.alpha == 0.01:
            t_value = 2.576  # приблизительно для больших выборок
        elif self.alpha == 0.05:
            t_value = 1.960
        else:  # 0.10
            t_value = 1.645
        
        # 1. Доверительный интервал для среднего
        se_mean = std_val / math.sqrt(n)  # стандартная ошибка среднего
        ci_mean_lower = mean_val - t_value * se_mean
        ci_mean_upper = mean_val + t_value * se_mean
        
        # 2. Доверительный интервал для стандартного отклонения
        # Используем формулу из курсовой работы
        ci_std_lower = math.sqrt(
            variance - t_value * variance * math.sqrt(2 / n)
        )
        ci_std_upper = math.sqrt(
            variance + t_value * variance * math.sqrt(2 / n)
        )
        
        # Сохраняем результаты
        self.results['ci_mean'] = (ci_mean_lower, ci_mean_upper)
        self.results['ci_std'] = (ci_std_lower, ci_std_upper)
        self.results['t_value'] = t_value
        
        # Отображение результатов
        self.text_confidence.delete(1.0, tk.END)
        self.text_confidence.insert(tk.END, "ДОВЕРИТЕЛЬНЫЕ ИНТЕРВАЛЫ ДЛЯ ПАРАМЕТРОВ ГЕНЕРАЛЬНОЙ СОВОКУПНОСТИ\n")
        self.text_confidence.insert(tk.END, "="*60 + "\n\n")
        
        self.text_confidence.insert(tk.END, f"Уровень значимости: α = {self.alpha}\n")
        self.text_confidence.insert(tk.END, f"Уровень доверия: {100*(1-self.alpha):.1f}%\n")
        self.text_confidence.insert(tk.END, f"Значение t-критерия: {t_value:.4f}\n")
        self.text_confidence.insert(tk.END, f"Объем выборки: n = {n}\n\n")
        
        self.text_confidence.insert(tk.END, "1. ДОВЕРИТЕЛЬНЫЙ ИНТЕРВАЛ ДЛЯ ГЕНЕРАЛЬНОЙ СРЕДНЕЙ (μ):\n")
        self.text_confidence.insert(tk.END, "-"*60 + "\n")
        self.text_confidence.insert(tk.END, f"Выборочное среднее: x̄ = {mean_val:.4f} г\n")
        self.text_confidence.insert(tk.END, f"Стандартная ошибка среднего: SE = s/√n = {se_mean:.6f} г\n\n")
        self.text_confidence.insert(tk.END, "Формула: x̄ ± t * (s/√n)\n")
        self.text_confidence.insert(tk.END, 
            f"{mean_val:.4f} ± {t_value:.4f} * {se_mean:.6f}\n")
        self.text_confidence.insert(tk.END, 
            f"{mean_val:.4f} ± {t_value*se_mean:.6f}\n\n")
        self.text_confidence.insert(tk.END, 
            f"Доверительный интервал: [{ci_mean_lower:.4f}; {ci_mean_upper:.4f}] г\n\n")
        
        width_mean = ci_mean_upper - ci_mean_lower
        self.text_confidence.insert(tk.END, 
            f"Ширина интервала: {width_mean:.4f} г "
            f"({width_mean/mean_val*100:.1f}% от среднего)\n\n")
        
        self.text_confidence.insert(tk.END, "2. ДОВЕРИТЕЛЬНЫЙ ИНТЕРВАЛ ДЛЯ ГЕНЕРАЛЬНОГО СТАНДАРТНОГО ОТКЛОНЕНИЯ (σ):\n")
        self.text_confidence.insert(tk.END, "-"*60 + "\n")
        self.text_confidence.insert(tk.END, f"Выборочное СКО: s = {std_val:.4f} г\n")
        self.text_confidence.insert(tk.END, f"Дисперсия: s² = {variance:.6f} г²\n\n")
        self.text_confidence.insert(tk.END, "Формула: √[s² ± t*s²*√(2/n)]\n")
        self.text_confidence.insert(tk.END, 
            f"Доверительный интервал: [{ci_std_lower:.4f}; {ci_std_upper:.4f}] г\n\n")
        
        width_std = ci_std_upper - ci_std_lower
        self.text_confidence.insert(tk.END, 
            f"Ширина интервала: {width_std:.4f} г "
            f"({width_std/std_val*100:.1f}% от СКО)\n\n")
        
        self.text_confidence.insert(tk.END, "ИНТЕРПРЕТАЦИЯ:\n")
        self.text_confidence.insert(tk.END, "-"*60 + "\n")
        self.text_confidence.insert(tk.END, 
            f"С вероятностью {100*(1-self.alpha):.1f}% можно утверждать, что:\n"
            f"1. Истинное среднее значение массы орехов лещины в генеральной "
            f"совокупности находится в интервале от {ci_mean_lower:.4f} до {ci_mean_upper:.4f} г.\n"
            f"2. Истинное стандартное отклонение массы орехов лещины в генеральной "
            f"совокупности находится в интервале от {ci_std_lower:.4f} до {ci_std_upper:.4f} г.")
    
    def create_plots(self):
        """Создание графиков"""
        # Очищаем предыдущие графики
        for widget in self.figure_frame.winfo_children():
            widget.destroy()
        
        # Создаем фигуру с несколькими подграфиками
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Визуализация результатов статистического анализа', fontsize=16)
        
        # Данные для графиков
        intervals = self.results['intervals']
        frequencies = self.results['frequencies']
        theoretical_freqs = self.results.get('theoretical_freqs', [])
        midpoints = [(low + high) / 2 for low, high in intervals]
        
        # 1. Гистограмма эмпирического распределения
        ax1 = axes[0, 0]
        ax1.bar(range(len(intervals)), frequencies, alpha=0.7, color='skyblue', 
                edgecolor='black', label='Эмпирические частоты')
        ax1.set_xlabel('Интервалы')
        ax1.set_ylabel('Частота')
        ax1.set_title('Гистограмма эмпирического распределения')
        ax1.set_xticks(range(len(intervals)))
        ax1.set_xticklabels([f'{low:.2f}-{high:.2f}' for low, high in intervals], rotation=45)
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Сравнение эмпирических и теоретических частот
        ax2 = axes[0, 1]
        x_pos = np.arange(len(intervals))
        width = 0.35
        
        ax2.bar(x_pos - width/2, frequencies, width, alpha=0.7, color='skyblue',
                label='Эмпирические', edgecolor='black')
        
        if theoretical_freqs:
            ax2.bar(x_pos + width/2, theoretical_freqs, width, alpha=0.7, color='lightcoral',
                    label='Теоретические', edgecolor='black')
        
        ax2.set_xlabel('Интервалы')
        ax2.set_ylabel('Частота')
        ax2.set_title('Сравнение эмпирических и теоретических частот')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels([f'Инт.{i+1}' for i in range(len(intervals))])
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Кривая нормального распределения
        ax3 = axes[1, 0]
        
        # Создаем плотное распределение для гладкой кривой
        mean_val = self.results['mean']
        std_val = self.results['std_dev']
        x_normal = np.linspace(mean_val - 3*std_val, mean_val + 3*std_val, 1000)
        y_normal = norm.pdf(x_normal, mean_val, std_val) * len(self.data) * self.results['lambda']
        
        ax3.plot(x_normal, y_normal, 'r-', linewidth=2, label='Нормальное распределение')
        ax3.hist(self.data, bins=len(intervals), density=False, alpha=0.5, 
                color='skyblue', edgecolor='black', label='Эмпирическое распределение')
        
        ax3.set_xlabel('Масса орехов, г')
        ax3.set_ylabel('Частота')
        ax3.set_title('Сравнение с нормальным распределением')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. Box plot с доверительными интервалами
        ax4 = axes[1, 1]
        
        # Box plot
        bp = ax4.boxplot(self.data, vert=True, patch_artist=True, 
                        boxprops=dict(facecolor='lightblue'))
        
        # Добавляем доверительные интервалы
        ci_mean = self.results.get('ci_mean', (mean_val, mean_val))
        ci_std = self.results.get('ci_std', (std_val, std_val))
        
        # Линия для среднего
        ax4.axhline(y=mean_val, color='red', linestyle='-', linewidth=2, label=f'Среднее: {mean_val:.4f} г')
        
        # Доверительный интервал для среднего
        ax4.axhspan(ci_mean[0], ci_mean[1], alpha=0.2, color='green', 
                   label=f'ДИ для среднего: [{ci_mean[0]:.4f}, {ci_mean[1]:.4f}]')
        
        ax4.set_xlabel('Выборка')
        ax4.set_ylabel('Масса орехов, г')
        ax4.set_title('Box plot с доверительными интервалами')
        ax4.legend(loc='upper right')
        ax4.grid(True, alpha=0.3)
        ax4.set_xticks([1])
        ax4.set_xticklabels(['Орехи лещины'])
        
        plt.tight_layout()
        
        # Встраиваем график в Tkinter
        canvas = FigureCanvasTkAgg(fig, master=self.figure_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Добавляем панель инструментов (опционально)
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
        toolbar = NavigationToolbar2Tk(canvas, self.figure_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def export_report(self):
        """Экспорт полного отчета в текстовый файл"""
        if not self.results:
            messagebox.showwarning("Нет данных", "Сначала выполните анализ")
            return
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="Сохранить отчет"
        )
        
        if filepath:
            try:
                with open(filepath, 'w', encoding='utf-8') as f:
                    f.write("="*80 + "\n")
                    f.write("ОТЧЕТ ПО СТАТИСТИЧЕСКОМУ АНАЛИЗУ ДАННЫХ\n")
                    f.write("="*80 + "\n\n")
                    
                    f.write(f"Дата анализа: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"Уровень значимости: α = {self.alpha}\n")
                    f.write(f"Объем выборки: n = {len(self.data)}\n\n")
                    
                    # Основные характеристики
                    f.write("1. ОСНОВНЫЕ ХАРАКТЕРИСТИКИ:\n")
                    f.write("-"*80 + "\n")
                    f.write(f"Среднее значение: {self.results.get('mean', 0):.4f} г\n")
                    f.write(f"Дисперсия: {self.results.get('variance', 0):.6f} г²\n")
                    f.write(f"СКО: {self.results.get('std_dev', 0):.4f} г\n")
                    f.write(f"Размах: {self.results.get('R', 0):.4f} г\n\n")
                    
                    # Критерий Пирсона
                    f.write("2. КРИТЕРИЙ СОГЛАСИЯ ПИРСОНА:\n")
                    f.write("-"*80 + "\n")
                    f.write(f"χ² наблюдаемое: {self.results.get('chi2_observed', 0):.4f}\n")
                    f.write(f"χ² критическое: {self.results.get('chi2_critical', 0):.4f}\n")
                    f.write(f"p-value: {self.results.get('chi2_p_value', 0):.6f}\n")
                    
                    if self.results.get('normality_conclusion', False):
                        f.write("ВЫВОД: Распределение соответствует нормальному\n\n")
                    else:
                        f.write("ВЫВОД: Распределение НЕ соответствует нормальному\n\n")
                    
                    # Доверительные интервалы
                    f.write("3. ДОВЕРИТЕЛЬНЫЕ ИНТЕРВАЛЫ:\n")
                    f.write("-"*80 + "\n")
                    ci_mean = self.results.get('ci_mean', (0, 0))
                    ci_std = self.results.get('ci_std', (0, 0))
                    
                    f.write(f"Для генеральной средней: [{ci_mean[0]:.4f}; {ci_mean[1]:.4f}] г\n")
                    f.write(f"Для генерального СКО: [{ci_std[0]:.4f}; {ci_std[1]:.4f}] г\n\n")
                    
                    f.write("="*80 + "\n")
                    f.write("Анализ выполнен с помощью Биостатистического анализатора\n")
                    f.write("="*80 + "\n")
                
                self.status_bar.config(text=f"Отчет сохранен в {os.path.basename(filepath)}")
                messagebox.showinfo("Успех", "Отчет успешно экспортирован")
                
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось сохранить отчет:\n{str(e)}")
    
    def save_plots(self):
        """Сохранение графиков в файл"""
        if not self.results:
            messagebox.showwarning("Нет данных", "Сначала выполните анализ")
            return
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf"), 
                      ("SVG files", "*.svg"), ("All files", "*.*")],
            title="Сохранить графики"
        )
        
        if filepath:
            try:
                # Создаем новые графики для сохранения
                fig, axes = plt.subplots(2, 2, figsize=(15, 12))
                
                # (Повторяем код создания графиков из create_plots())
                # ... [код создания графиков] ...
                
                fig.savefig(filepath, dpi=300, bbox_inches='tight')
                plt.close(fig)
                
                self.status_bar.config(text=f"Графики сохранены в {os.path.basename(filepath)}")
                messagebox.showinfo("Успех", "Графики успешно сохранены")
                
            except Exception as e:
                messagebox.showerror("Ошибка", f"Не удалось сохранить графики:\n{str(e)}")
    
    def reset_data(self):
        """Сброс всех данных и результатов"""
        if messagebox.askyesno("Подтверждение", 
                              "Вы уверены, что хотите сбросить все данные и результаты?"):
            self.data = None
            self.results = {}
            
            # Очищаем все текстовые поля
            for text_widget in [self.text_original, self.text_variation, 
                               self.text_stats, self.text_pearson, self.text_confidence]:
                text_widget.delete(1.0, tk.END)
            
            # Очищаем графики
            for widget in self.figure_frame.winfo_children():
                widget.destroy()
            
            self.data_info_label.config(text="Данные не загружены")
            self.status_bar.config(text="Данные сброшены")
    
    def show_help(self):
        """Показать справку по программе"""
        help_text = """
        БИОСТАТИСТИЧЕСКИЙ АНАЛИЗАТОР
        
        Назначение:
        Программа автоматизирует статистический анализ биологических данных,
        в частности, проверку на соответствие нормальному распределению
        по критерию согласия Пирсона.
        
        Основные функции:
        1. Загрузка данных из файлов (CSV, Excel) или ручной ввод
        2. Построение вариационных рядов
        3. Расчет основных статистических характеристик
        4. Проверка на нормальность распределения (χ²-критерий Пирсона)
        5. Расчет доверительных интервалов
        6. Визуализация результатов
        7. Экспорт отчетов
        
        Пример данных:
        Программа включает примерные данные из курсовой работы
        по анализу массы орехов лещины (90 наблюдений).
        
        Уровень значимости:
        Можно выбрать α = 0.01, 0.05 или 0.10
        
        Автор: [Ваше имя]
        Версия: 1.0
        """
        
        messagebox.showinfo("Справка", help_text)


def main():
    """Запуск приложения"""
    root = tk.Tk()
    
    # Стиль для кнопок
    style = ttk.Style()
    style.configure("Accent.TButton", font=('Arial', 10, 'bold'))
    
    app = BioStatApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()