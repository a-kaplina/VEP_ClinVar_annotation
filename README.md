# VEP_ClinVar_annotation

Скрипт  обрабатывает VCF файл **(в сборке hg38)**, добавляет к нему аннотации (VEP, ClinVar, PharmGKB), фильтрует по клинической значимости и создаёт итоговые таблицы для анализа.
В папке со скриптом должны находиться необходимые файлы - make_long_table_script.py и pharmgkb_map_clean_CAT_all.txt

Применение (нужно указать путь к папке, в которой находится VCF файл для аннотации **в сборке hg38**): \
bash VEP_Clinvar_annot_script.sh ~/absolute/path/VCF

**Для работы скрипта необходимо:** 

**Предварительно установить VEP**
(https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html):
```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
```
```bash
git pull
git checkout release/115
perl INSTALL.pl
```

 **Проверить, есть ли файлы сборки hg38 в папке VEP:**
```bash
ls ~/.vep/homo_sapiens
```

Внутри должна быть папка: 115_GRCh38

- Если ее нет, сначала надо скачать homo_sapiens_vep_115_GRCh38.tar.gz для VEP:
```bash
cd $HOME/.vep/homo_sapiens
curl -O https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
tar xzf homo_sapiens_vep_115_GRCh38.tar.gz
mv ~/.vep/homo_sapiens/homo_sapiens/115_GRCh38 ~/.vep/homo_sapiens/
```
Должно быть так:
```bash
ls $HOME/.vep/homo_sapiens
```
115_GRCh37  115_GRCh38  homo_sapiens_vep_115_GRCh38.tar.gz

**Далее нужно скачать primary assembly:**
```bash
cd HOME/.vep/homo_sapiens/115_GRCh38
wget ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
- Надо перезаписать ее в bgz
```bash
zcat ~/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | bgzip -c > ~/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bgz
```

**Скачать ClinVar для hg 38:**
```bash
mkdir ~/clinvar38
cd ~/clinvar38
 wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
bcftools index clinvar.vcf.gz
```
**Схематичное отображение этапов работы**
<br>
<img width="1181" height="902" alt="VEP_diagram drawio" src="https://github.com/user-attachments/assets/03922237-236e-491c-aaac-c72597cbbfef" />

**Описание работы скрипта:** \
Проверка входных данных

    Проверяет, указана ли папка с исходными файлами

    Создаёт внутри неё структуру папок: input/ (временная папка) и output/

Шаги 1 и 2: Сортировка и нормализация VCF файла

    Сортирует исходный VCF файл, разделяет мультиаллельные варианты на отдельные записи, приводит варианты к стандартизированному виду

    Создаёт файл: *_norm.vcf.gz

Шаг 3: Cоздание mapping файла для удаления chr из названий хромосом в VCF

Шаг 4: Удаление префикса "chr"

    Убирает приставку chr из названий хромосом (chr1 → 1, chrX → X) для совместимости с базой ClinVar

    Создаёт файл: *_nochr.vcf.gz

Шаг 5: Очистка промежуточных ( _norm) файлов 

Шаг 6: VEP аннотация + ClinVar

    Запускает VEP (Ensembl Variant Effect Predictor) для добавления:

        Функциональных характеристик вариантов (missense, nonsense и т.д.)

        Предсказаний SIFT и PolyPhen 

        Частот из 1000 Genomes и gnomAD 

        HGVS номенклатуры

        Информации о генах, транскриптах, доменах белков

    Затем добавляет данные из ClinVar (клиническая значимость, фенотипы, статус рецензирования)

    Создаёт файл: *_annotated.vcf.gz

Шаг 7: Удаление временного _nochr файла

Шаг 8: Добавление PharmGKB аннотаций с помощью предварительно подготовленного файла pharmgkb_map_clean_CAT_all.txt (создан из ClinicalVariants.tsv с сайта ClinPGX с добавлением вариантов, полученных при помощи аннотации PharmCAT - указаны, как PharmCAT)

    Добавляет информацию из базы PharmGKB (уровни доказательств 1A, 1B, 2A, 2B, 3, 4, а также варианты, полученные ранее при аннотации через PharmCAT)

    Создаёт файл: *_annotated_Pharm.vcf.gz

    Подсчитывает количество вариантов из PharmGKB с различными уровнями (1-4, PharmCAT)

Шаг 9: Удаление промежуточного аннотированного файла

    Удаляет *_annotated.vcf.gz, оставляя только версию с PharmGKB

Шаг 10: Фильтрация по генотипу

    Убирает варианты, для которых все образцы гомозиготны по референсному аллелю

    Создаёт файл: *_real_annotated_Pharm.vcf.gz

Шаг 11: Статистика по патогенности

    Создаёт файлы _clnsig_counts.txt с подсчётом различных классов клинической значимости и review status (для файла до фильтрации и после удаления гомозигот по референсу)

Шаг 12: Выделение патогенных вариантов

    Оставляет только варианты с метками: Pathogenic, Likely_pathogenic и их комбинации

    Создаёт файл: *__all_pathogenic.vcf.gz

Шаг 13: Фильтрация по качеству рецензирования (патогенные)

    Оставляет только патогенные варианты с хорошим статусом рецензирования:

        2 звезды (множество подтверждений, нет конфликтов)

        3 звезды (экспертная панель)

        4 звезды (клинические рекомендации)

        1 звезда (хотя бы одно подтверждение с критериями)

    Создаёт файл: *_good_q_pathogenic.vcf.gz

Шаг 14: Выделение VUS вариантов

    Фильтрует аннотированный файл отдельно и оставляет только варианты с Uncertain Significance (неопределённая значимость)

    Создаёт файл: *_all_vus.vcf.gz

Шаг 15: Фильтрация VUS по качеству рецензирования

    Оставляет только VUS варианты с хорошим статусом рецензирования

    Создаёт файл: *_good_q_vus.vcf.gz

Шаг 16 и 17: Создание таблиц при помощи скрипта make_long_table_script.py

    Запускает Python скрипт для конвертации VCF в табличный формат

    Создаёт таблицы:

        *_Result_pat_fin_annotation.tsv - для патогенных вариантов

        *_Result_vus_fin_annotation.tsv - для VUS вариантов

    В таблицах содержатся все поля: частоты, предсказания, фенотипы, лекарства и т.д.
